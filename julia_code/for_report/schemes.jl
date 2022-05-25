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
    B, D, a, b, x = (zeros(Nx) for i in 1:5)
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
        tridiagonal_algorithm!(a, b, x, A, B, C, D) # Метод прогонки
        u[:, j + 1] .= x
    end

    return u
end