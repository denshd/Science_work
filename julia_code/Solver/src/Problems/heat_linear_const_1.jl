export solve_PDE

function solve_PDE(sys::Heat_linear_const_1, g::Grid, method::Symbol; σ=0.5)
    # Матрица, хранящая решение разностной задачи
    u = zeros(g.Nx, g.Nt)

    # Учёт начальных условий (действие оператора P_h)
    initial_cond!(u, sys, g)

    # Учёт граничных условий
    boundary_cond_1!(u, sys, g)

    # Решение уравнения в зависимости от выбранного метода
    if (method == :explicit)
        explicit!(u, sys, g)
    elseif (method == :σ_implicit)
        σ_implicit!(u, sys, g, σ)
    end

    return Solution(u, g)
end


"""
Явная схема на равномерной сетке
"""
function explicit!(u::Matrix{Float64}, sys::Heat_linear_const_1, g::UniformGrid)
    c = g.τ / (g.h^2)
    # Проверка схемы на устойчивость
    if (c > 0.5)
        throw(ErrorException("Явная схема неустойчива!, k=$c"))
    end

    # Аппроксимация неоднородности уравнения
    ϕ = zeros(g.Nx, g.Nt)
    for i = 1:g.Nx
        for j = 1:g.Nt
            ϕ[i, j] = sys.f(g.x[i], g.t[j])
        end
    end

    # Основной цикл явной схемы
    for j = 1:g.Nt - 1
        for i = 2:g.Nx - 1
            u[i, j+1] = u[i, j] + c * (u[i+1, j] - 2 * u[i, j] + u[i-1, j]) + ϕ[i, j]
        end
    end
end


"""
Явная схема на неравномерной сетке
"""
function explicit!(u::Matrix{Float64}, sys::Heat_linear_const_1, g::GeneralGrid)
    # c = max(g.τ...) / (min(g.h...)^2)
    # Проверка схемы на устойчивость
    if (c > 0.5)
        throw(ErrorException("Явная схема неустойчива!, k=$c"))
    end

    # Аппроксимация неоднородности уравнения
    ϕ = zeros(g.Nx, g.Nt)
    for i = 1:g.Nx
        for j = 1:g.Nt
            ϕ[i, j] = sys.f(g.x[i], g.t[j])
        end
    end

    # Основной цикл явной схемы
    for j = 1:g.Nt - 1
        for i = 2:g.Nx - 1
            u[i, j+1] = u[i, j] + (2 * g.τ[j] / (g.h[i] + g.h[i-1])) * (
                (u[i+1, j] - u[i, j]) / g.h[i] -
                (u[i, j] - u[i-1, j]) / g.h[i-1]
            )
        end
    end
end


"""
Неявная σ-схема на равномерной сетке
"""
function σ_implicit!(u::Matrix{Float64}, sys::Heat_linear_const_1, g::UniformGrid, σ::Real)
    k = 0.5 - g.h^2/(4*g.τ)
    # Проверка на устойчивость
    if (σ < 0.5)
        throw(ErrorException("Неявная схема с параметром σ=$σ неустойчива!, k=$k"))
    end

    # Аппроксимация неоднородности уравнения
    ϕ = zeros(g.Nx, g.Nt)
    for i = 1:g.Nx
        for j = 1:g.Nt
            ϕ[i, j] = sys.f(g.x[i], g.t[j])
        end
    end

    # Основной цикл неявной схемы

    ## нижняя диагональ
    A = -σ * g.τ * ones(g.Nx - 1)
    A[end] = 0
    
    ## главная диагональ
    B = (g.h^2 + 2*σ*g.τ) * ones(g.Nx)
    B[1] = B[end] = 1

    # верхняя диагональ
    C = -σ * g.τ * ones(g.Nx - 1)
    C[1] = 0

    D = zeros(g.Nx) # столбец свободных членов

    # Вспомогательные векторы для прогонки
    α = zeros(g.Nx)
    β = zeros(g.Nx)
    temp_x = zeros(g.Nx)

    for j = 1 : g.Nt-1
        for i = 2 : g.Nx-1
            D[i] = g.h^2 * u[i, j] + g.τ*(1-σ) * (
                u[i-1, j] - 2 * u[i, j] + u[i+1, j]
            ) + g.τ * g.h^2 * ϕ[i, j+1]
        end
        D[1] = u[1, j+1]
        D[end] = u[end, j+1]

        tridiagonal_matrix_algorithm!(α, β, temp_x, A, B, C, D)
        u[:, j+1] .= temp_x

        # метод прогонки
        # @time u[:, j+1] = M\D
    end
end


"""
Неявная σ-схема на неравномерной сетке
"""
function σ_implicit!(u::Matrix{Float64}, sys::Heat_linear_const_1, g::GeneralGrid, σ::Real)
    k = 0.5 - max(g.h...)^2/(4*min(g.τ...))
    # Проверка на устойчивость
    if (σ < 0.5)
        throw(ErrorException("Неявная схема с параметром σ=$σ неустойчива!, k=$k"))
    end

    # Аппроксимация неоднородности уравнения
    ϕ = zeros(g.Nx, g.Nt)
    for i = 1:g.Nx
        for j = 1:g.Nt
            ϕ[i, j] = sys.f(g.x[i], g.t[j])
        end
    end

    # Основной цикл неявной схемы

    ## нижняя диагональ
    A::Vector{Float64} = vcat([-2*σ / (g.h[i] * (g.h[i+1] + g.h[i])) for i in 1:(g.Nx - 2)], [0])

    # Главная диагональ
    B = zeros(g.Nx)
    B[1] = 1
    B[end] = 1

    # верняя диагональ
    C::Vector{Float64} = vcat([0], [-2 * σ / (g.h[i] * (g.h[i] + g.h[i - 1])) for i in 2:(g.Nx - 1)])

    D = zeros(g.Nx) # столбец свободных членов

    # Вспомогательные векторы для прогонки
    α = zeros(g.Nx)
    β = zeros(g.Nx)
    temp_x = zeros(g.Nx)

    for j = 1 : g.Nt-1
        for i = 2 : g.Nx-1
            D[i] = u[i, j] / g.τ[j] +
            (1 - σ) * (
                (u[i + 1, j] - u[i, j]) / g.h[i] -
                (u[i, j] - u[i - 1, j]) / g.h[i - 1]
            ) *
            2 / (g.h[i] + g.h[i - 1])

            B[i] = 1 / g.τ[j] +
            2 * σ / (g.h[i] * (g.h[i] + g.h[i - 1])) +
            2 * σ / (g.h[i - 1] * (g.h[i] + g.h[i - 1]))
        end
        D[1] = u[1, j+1]
        D[end] = u[end, j+1]

        tridiagonal_matrix_algorithm!(α, β, temp_x, A, B, C, D)
        u[:, j+1] .= temp_x

        # метод прогонки
        # @time u[:, j+1] = M\D
    end
end