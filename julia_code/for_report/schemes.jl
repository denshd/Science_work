function explicit_scheme(f, u0, mu1, mu2, T, Nx, Nt)
    x = range(0, 1, length = Nx) # Точки секти по оси `x`
    t = range(0, T, length = Nt) # Точки сетки по оси `t`
    c = step(t) / step(x)^2 # Коэффициент Куранта

    if (c > 0.5)
        throw(ErrorException("Явная схема неустойчива, c = $(c) > 0.5"))
    end

    u = zeros(Nx, Nt) # Матрица, хранящая решение разностной задачи

    u[:, 1] .= u0.(x) # Учёт начальных условий
    u[1, :] .= mu1.(t) # Учёт граничных условий на левом конце
    u[end, :] .= mu2.(t) # Учёт граничных условий на правом конце
 
    for j in 1:(Nt - 1) # Цикл по всем слоям
        for i in 2:(Nx - 1) # Явно выражаем значения функции на новом слое через значения на старом
            u[i, j + 1] = u[i, j] + c * (u[i + 1, j] - 2 * u[i, j] + u[i - 1, j]) + f(x[i], t[j])
        end
    end

    return u
end


# Оптимальный по памяти алгоритм прогонки в случае, если его нужно запускать много раз
function tridiagonal_algorithm!(
    α::Vector{<:T}, β::Vector{<:T}, x::Vector{<:T},
    A::Vector{<:V}, B::Vector{<:V}, C::Vector{<:V}, D::Vector{<:V}
    ) where {T <: Number, V <: Number}

    α[1] = - C[1] / B[1]
    β[1] = D[1] / B[1]
    @inbounds for i=2:(length(B)-1)
        α[i] = - C[i] / (B[i] + A[i-1] * α[i-1])
        β[i] = (D[i] - A[i-1] * β[i-1]) / (B[i] + A[i-1] * α[i-1])
    end
    β[end] = (D[end] - A[end] * β[end-1]) / (B[end] + A[end] * α[end-1])

    # обратный ход
    x[end] = β[end]
    for i = length(B)-1:-1:1
        x[i] = α[i]*x[i+1] + β[i]
    end
end


function implicit_scheme(k, f, u0, mu1, mu2, T, Nx, Nt)
    x = range(0, 1, length = Nx) # Точки секти по оси `x`
    dx = step(x)
    t = range(0, T, length = Nt) # Точки сетки по оси `t`
    dt = step(t)

    u = zeros(Nx, Nt) # Матрица, хранящая решение разностной задачи
    u[:, 1] .= u0.(x) # Учёт начальных условий
    u[1, :] .= mu1.(t) # Учёт граничных условий на левом конце
    u[end, :] .= mu2.(t) # Учёт граничных условий на правом конце

    # Временные переменные для алгоритма прогонки
    A, C = (zeros(Nx) for i in 1:2)
    B, D, a, b, xt = (zeros(Nx) for i in 1:5)
    B[1] = B[end] = 1
    C[1] = A[end] = 0

    for j in 1:(Nt - 1) # Цикл по временным слоям
        # Обновление коэффициентов тридиагональной матрицы
        D[1] = u[1, j + 1]
        D[end] = u[end, j + 1]
        for i = 2:(Nx - 1)
            A[i - 1] = -dt * k(0.5 * (u[i - 1, j] + u[i, j]))
            B[i] = dx^2 + dt * (
                k(0.5 * (u[i, j] + u[i + 1, j])) +
                k(0.5 * (u[i - 1, j] + u[i, j]))
            )
            C[i] = -dt * k(0.5 * (u[i, j] + u[i + 1, j]))
            D[i] = dx^2 * u[i, j] + dt * dx^2 * f(x[i], t[j])
        end
        tridiagonal_algorithm!(a, b, xt, A, B, C, D) # Метод прогонки
        u[:, j + 1] .= xt
    end

    return u
end


function local_scheme(
    k1, k2, # коэффициенты теплопроводности
    f, # неоднородность уравнения
    u0, # начальные данные $u_0(x, y)$
    mu_1, mu1, mu_2, mu2, # граничные услвоия
    Lx, Ly, T, # параметры цилиндра
    Nx, Ny, Nt # число точек сетки 
    )

    # Создание сетки
    x = range(0, Lx, length = Nx)
    dx = step(x)
    y = range(0, Ly, length = Ny)
    dy = step(y)
    t = range(0, T, length = Nt)
    dt = step(t)

    u = zeros(Nx, Ny, Nt) # Массив, хранящий решение
    
    # Учёт начальных условий
    for i in 1:Nx, j in 1:Ny u[i, j, 1] = u0(x[i], y[j]) end

    # Временные переменные для прогонки
    A, C = (Tuple(zeros(N - 1) for N in [Nx Ny]) for i in 1:2)
    B, D, a, b, xt = (Tuple(zeros(N) for N in [Nx Ny]) for i in 1:5)
    B[1][1] = B[1][end] = B[2][1] = B[2][end] = 1
    C[1][1] = C[2][1] = A[1][end] = A[2][end] = 0
    w = zeros(Nx, Ny)


    @inbounds for j in 1:(Nt - 1) # Цикл по временным слоям
        @inbounds for k in 1:Ny # Сначала решаем задачу вдоль оси $x$ для всех $y_k$
            @inbounds for i in 2:(Nx - 1) # Создание трёхдиагональной матрицы
                A[1][i - 1] = -k1(0.5 * (u[i - 1, k, j] + u[i, k, j])) / dx^2
                B[1][i] = 1/dt + (
                    k1(0.5 * (u[i - 1, k, j] + u[i, k, j])) +
                    k1(0.5 * (u[i + 1, k, j] + u[i, k, j]))
                ) / dx^2
                C[1][i] = -k1(0.5 * (u[i + 1, k, j] + u[i, k, j])) / dx^2
                D[1][i] = u[i, k, j] / dt + 0.5 * f(x[i], y[k], t[j] + 0.5 * dt)
            end
            D[1][1] = mu_1(y[k], t[j] + 0.5 * dt)
            D[1][end] = mu1(y[k], t[j] + 0.5 * dt)

            tridiagonal_algorithm!(a[1], b[1], xt[1], A[1], B[1], C[1], D[1]) # Метод прогонки
            w[:, k] .= xt[1]
        end
        @inbounds for i in 1:Nx # Теперь решаем задачу вдоль оси $y$ для всех $x_i$
            @inbounds for k in 2:(Ny - 1) # Вся суть далее та же самая
                A[2][k - 1] = -k2(0.5 * (w[i, k] + w[i, k - 1])) / dy^2
                B[2][k] = 1/dt + (
                    k2(0.5 * (w[i, k] + w[i, k - 1])) +
                    k2(0.5 * (w[i, k + 1] + w[i, k]))
                ) / dy^2
                C[2][k] = -k2(0.5 * (w[i, k] + w[i, k + 1])) / dy^2
                D[2][k] = w[i, k] / dt + 0.5 * f(x[i], y[k], t[j + 1])
            end
            D[2][1] = mu_2(x[i], t[j + 1])
            D[2][end] = mu2(x[i], t[j + 1])
            tridiagonal_algorithm!(a[2], b[2], xt[2], A[2], B[2], C[2], D[2]) # Метод прогонки
            u[i, :, j + 1] .= xt[2]
        end
    end

    return u
end


# function local_scheme(
#     k1, k2, # коэффициенты теплопроводности
#     f, # неоднородность уравнения
#     u0, # начальные данные $u_0(x, y)$
#     mu_1, mu1, mu_2, mu2, # граничные услвоия
#     Lx, Ly, T, # параметры цилиндра
#     Nx, Ny, Nt # число точек сетки 
#     )

#         # Создание сетки
#     x = range(0, Lx, length = Nx)
#     dx = step(x)
#     y = range(0, Ly, length = Ny)
#     dy = step(y)
#     t = range(0, T, length = Nt)
#     dt = step(t)

#     u = zeros(Nx, Ny, Nt)

#     # Начальные данные
#     for i in 1:Nx, j in 1:Ny u[i, j, 1] = u0(x[i], y[j]) end

#     # Временное решение (см. алгоритм)
#     w = zeros(Nx, Ny)

#     # Диагонали в одномерных задачах
#     A::Tuple{Vector{Float64}, Vector{Float64}} = (
#         Vector{Float64}(undef, Nx - 1),
#         Vector{Float64}(undef, Ny - 1)
#     )
#     B::Tuple{Vector{Float64}, Vector{Float64}} = (
#         Vector{Float64}(undef, Nx),
#         Vector{Float64}(undef, Ny)
#     )
#     C::Tuple{Vector{Float64}, Vector{Float64}} = (
#         Vector{Float64}(undef, Nx - 1),
#         Vector{Float64}(undef, Ny - 1)
#     )
#     D::Tuple{Vector{Float64}, Vector{Float64}} = (
#         Vector{Float64}(undef, Nx),
#         Vector{Float64}(undef, Ny)
#     )
#     B[1][1] = B[1][end] = B[2][1] = B[2][end] = 1
#     C[1][1] = C[2][1] = A[1][end] = A[2][end] = 0

#     # Временные векторы для прогонки
#     α::Tuple{Vector{Float64}, Vector{Float64}} = (
#         Vector{Float64}(undef, Nx),
#         Vector{Float64}(undef, Ny)
#     )
#     β::Tuple{Vector{Float64}, Vector{Float64}} = (
#         Vector{Float64}(undef, Nx),
#         Vector{Float64}(undef, Ny)
#     )
#     temp_x::Tuple{Vector{Float64}, Vector{Float64}} = (
#         Vector{Float64}(undef, Nx),
#         Vector{Float64}(undef, Ny)
#     )

#     # Цикл по временным слоям
#     for j = 1 : Nt - 1

#         t_half = t[j] + dt / 2

#         # Сначала решаем задачу вдоль оси x1 (для всех i2)
#         for i2 = 1 : Ny

#             # Создаём трёхдиагональную матрицу
#             for i1 = 2 : Nx - 1
#                 A[1][i1 - 1] = -(1 /dx^2) * k1(0.5 * (u[i1 - 1, i2, j] + u[i1, i2, j]))
#                 B[1][i1] = 1 / dt + (1 /dx^2) * (
#                     k1(0.5 * (u[i1 - 1, i2, j] + u[i1, i2, j])) +
#                     k1(0.5 * (u[i1 + 1, i2, j] + u[i1, i2, j]))
#                 )
#                 C[1][i1] = -(1 /dx^2) * k1(0.5 * (u[i1 + 1, i2, j] + u[i1, i2, j]))
#                 D[1][i1] = (1 / dt) * u[i1, i2, j]
#             end
#             B[1][1] = B[1][end] = 1
#             C[1][1] = A[1][end] = 0
#             D[1][1] = mu_1(y[i2], t_half)
#             D[1][end] = mu1(y[i2], t_half)

#             # Получаем временный вектор v_(1)
#             tridiagonal_algorithm!(α[1], β[1], temp_x[1], A[1], B[1], C[1], D[1])
#             w[:, i2] .= temp_x[1]
#         end

#         # Теперь решаем локально-одномерную задачу вдоль оси x2 (для всех i1)
#         for i1 = 1 : Nx

#             # Создаём трёхдиагональную матрицу
#             for i2 = 2 : Ny - 1
#                 A[2][i2 - 1] = -(1 / dy^2) * k2(0.5 * (w[i1, i2] + w[i1, i2 - 1]))
#                 B[2][i2] = 1 / dt + (1 / dy^2) * (
#                     k2(0.5 * (w[i1, i2] + w[i1, i2 - 1])) +
#                     k2(0.5 * (w[i1, i2 + 1] + w[i1, i2]))
#                 )
#                 C[2][i2] = -(1 / dy^2) * k2(0.5 * (w[i1, i2] + w[i1, i2 + 1]))
#                 D[2][i2] = (1 / dt) * w[i1, i2]
#             end
#             B[2][1] = B[2][end] = 1
#             C[2][1] = A[2][end] = 0
#             D[2][1] = mu_2(x[i1], t[j+1])
#             D[2][end] = mu2(x[i1], t[j+1])
#             # D[2][1] = u[i1, 1, j + 1]
#             # D[2][end] = u[i1, end, j + 1]

#             # Получаем решение на временном срезе j + 1 для x[1] = x[1][i]
#             tridiagonal_algorithm!(α[2], β[2], temp_x[2], A[2], B[2], C[2], D[2])
#             u[i1, :, j + 1] .= temp_x[2]
#         end
#     end

#     return u
# end