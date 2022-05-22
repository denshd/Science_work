export solve_PDE

function solve_PDE(sys::Heat_quazilinear_1, g::UniformGrid)
    # Представление решения
    u = zeros(g.Nx, g.Nt)

    # Начальные/граничные условия
    initial_cond!(u, sys, g)
    boundary_cond_1!(u, sys, g)

    # Векторы для прогонки
    A = zeros(g.Nx - 1)
    B = zeros(g.Nx)
    C = zeros(g.Nx - 1)
    D = zeros(g.Nx)

    # Вспомогательные векторы для прогонки
    α = zeros(g.Nx)
    β = zeros(g.Nx)
    temp_x = zeros(g.Nx)

    # цикл по временным слоям
    for j = 1:g.Nt-1

        D[1] = u[1, j+1]
        B[1] = 1
        C[1] = 0
        for i = 2:g.Nx-1
            A[i-1] = -g.τ * sys.k((u[i-1, j] + u[i, j]) / 2)
            B[i] = g.h^2 + g.τ * sys.k((u[i, j] + u[i+1, j]) / 2) + g.τ * sys.k((u[i-1, j] + u[i, j]) / 2)
            C[i] = - g.τ * sys.k((u[i, j] + u[i+1, j]) / 2)
            D[i] = g.τ * g.h^2 * sys.f(g.x[i], g.t[j], u[i, j]) + g.h^2 * u[i, j]
        end
        D[end] = u[end, j+1]
        B[end] = 1
        A[end] = 0
        
        tridiagonal_matrix_algorithm!(α, β, temp_x, A, B, C, D)
        u[:, j+1] .= temp_x
    end

    return Solution(u, g);
end


function solve_PDE(sys::Heat_quazilinear_1, g::GeneralGrid)
    # Представление решения
    u = zeros(g.Nx, g.Nt)

    # Начальные/граничные условия
    initial_cond!(u, sys, g)
    boundary_cond_1!(u, sys, g)

    # Векторы для прогонки
    A = zeros(g.Nx - 1)
    B = zeros(g.Nx)
    C = zeros(g.Nx - 1)
    D = zeros(g.Nx)

    # Вспомогательные векторы для прогонки
    α = zeros(g.Nx)
    β = zeros(g.Nx)
    temp_x = zeros(g.Nx)

    # "Средний шаг"
    ħ = 0.5 .* (g.h[1:end-1] .+ g.h[2:end])

    # Цикл по временным слоям
    for j = 1:g.Nt-1
        D[1] = u[1, j+1]
        B[1] = 1
        C[1] = 0
        for i = 2:g.Nx-1
            A[i-1] = - sys.k((u[i-1, j] + u[i, j]) / 2) / (ħ[i-1] * g.h[i-1])
            B[i] = 1 / g.τ[j] + sys.k((u[i+1, j] + u[i, j]) / 2) / (ħ[i-1] * g.h[i]) + sys.k((u[i-1, j] + u[i, j]) / 2) / (ħ[i-1] * g.h[i-1]) 
            C[i] = - sys.k((u[i, j] + u[i+1, j]) / 2) / (ħ[i-1] * g.h[i])
            D[i] = sys.f(g.x[i], g.t[j], u[i, j]) + u[i, j] / g.τ[j]
        end
        D[end] = u[end, j+1]
        B[end] = 1
        A[end] = 0
        
        tridiagonal_matrix_algorithm!(α, β, temp_x, A, B, C, D)
        u[:, j+1] .= temp_x
    end

    return Solution(u, g);
end


# function solve_PDE(sys::Heat_quazilinear_1, g::AGrid, subdivision_level=2)
#     u = zeros(g.Nx, g.Nt)

#     initial_cond!(u, sys, g)
#     boundary_cond_1!(u, sys, g)

#     # Основной цикл неявной схемы

#     # нижняя диагональ
#     A = zeros(g.Nx - 1)
#     B = zeros(g.Nx)
#     C = zeros(g.Nx - 1)
#     D = zeros(g.Nx)

#     # Вспомогательные векторы для прогонки
#     α = zeros(g.Nx)
#     β = zeros(g.Nx)
#     temp_x = zeros(g.Nx)


#     # Цикл по слоям
#     for j = 1:g.Nt-1

#         # решаем один слой на coarse сетке
#         solve_one_layer!(u, j, sys, g.coarse, A, B, C, D, α, β, temp_x)

#         # Запускаем алгоритм дробёжки
#         regrid!(u, j, sys, g, 0)

#         # ищем, где дробить сетку
#     end

#     return Solution(u, g);
# end


# function regrid!(
#         u::Matrix{Float64},
#         j::Int64,
#         sys::Heat_quazilinear_1,
#         g::AGrid, current_level::Int64, super_parameter::Int64=4
#     )

#     if (current_level == g.subdivision_level)
#         return;
#     end

#     # Находим максимум градиента
#     i0 = argmax(u[3:end] .- u[1:end-2])

#     # Определяем границы fine-сетки
#     temp_Nx = min(g.Nx - i0, i0, div(g.coarse.Nx, super_parameter))
#     new_Nx = 2 + 4 * temp_Nx - 1

#     # Создаём новую Uniform-сетку
#     g_new = UniformGrid((g.coarse.x[i0-temp_Nx], g.coarse.x[i0+temp_Nx]), (0, 0), new_Nx, 0)

#     # Решаем один слой на новой Uniform-сетке







# function solve_one_layer!(
#         u::Matrix{Float64},
#         j::Int64,
#         sys::Heat_quazilinear_1,
#         g::UniformGrid,
#         A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64}, D::Vector{Float64},
#         α::Vector{Float64}, β::Vector{Float64}, temp_x::Vector{Float64}
#     )

#     D[1] = u[1, j+1]
#     B[1] = 1
#     C[1] = 0
#     for i = 2:g.Nx-1
#         A[i-1] = -g.τ * sys.k((u[i-1, j] + u[i, j]) / 2)
#         B[i] = g.h^2 + g.τ * sys.k((u[i, j] + u[i+1, j]) / 2) + g.τ * sys.k((u[i-1, j] + u[i, j]) / 2)
#         C[i] = - g.τ * sys.k((u[i, j] + u[i+1, j]) / 2)
#         D[i] = g.τ * g.h^2 * sys.f(g.x[i], g.t[j], u[i, j]) + g.h^2 * u[i, j]
#     end
#     D[end] = u[end, j+1]
#     B[end] = 1
#     A[end] = 0
    
#     tridiagonal_matrix_algorithm!(α, β, temp_x, A, B, C, D)
#     u[:, j+1] .= temp_x
# end
