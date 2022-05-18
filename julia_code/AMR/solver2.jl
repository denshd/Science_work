include("solver.jl")

"""
Равномерная двумерная сетка из клеток и граней 
(Представляется просто как декартово произведение одномерных)
"""
struct UniformGrid2{T <: AbstractRange{Cdouble}}
    g::Tuple{UniformGrid{T}, UniformGrid{T}}
    

    """
    Создаёт двумерную равномерную сетку на `[0, 1] × [0, 1]` с числом клеток `Ncells::Tuple{Integer, Integer}`
    """
    function UniformGrid2(Ncells::Tuple{Integer, Integer})
        g = Tuple(UniformGrid(Ncells[i]) for i = 1:2)
        new{typeof(g[1].cells)}(g)
    end


    """
    Создаёт двумерную равномерную сетку на `[0, 1] × [0, 1]` с шагами `h::Tuple{Real, Real}`
    """
    function UniformGrid2(h::Tuple{Real, Real})
        g = Tuple(UniformGrid(h[i]) for i = 1:2)
        new{typeof(g[1].cells)}(g)
    end


    """
    Сооздаёт двумерную равномерную сетку на `[0, L[1]] × [0, L[2]]` с числом точек `Ncells::Tuple{Integer, Integer}`
    """
    function UniformGrid2(L::Tuple{Real, Real}, Ncells::Tuple{Integer, Integer})
        g = Tuple(UniformGrid(L[i], Ncells[i]) for i = 1:2)
        new{typeof(g[1].cells)}(g)
    end


    """
    Создаёт двумерную равномерную сетку на `[0, L[1]] × [0, L[2]]` с шагами `h::Tuple{Real, Real}`
    """
    function UniformGrid2(L::Tuple{Real, Real}, h::Tuple{Real, Real})
        g = Tuple(UniformGrid(L[i], h[i]) for i = 1:2)
        new{typeof(g[1].cells)}(g)
    end


    """
    Создаёт двумерную равномерную сетку на `[x_b[1][1], x_b[1][2]] × [x_b[2][1], x_b[2][2]` с шагмаи `Nsteps::Tuple{Integer, Integer}`
    """
    function UniformGrid2(x_b::Tuple{Tuple{Real, Real}, Tuple{Real, Real}}, Ncells::Tuple{Integer, Integer})
        g = Tuple(UniformGrid(x_b[i], Ncells[i]) for i = 1:2)
        new{typeof(g[1].cells)}(g)
    end
end


struct UniformGrid2Solution{T <: AbstractRange{Cdouble}}
    spacial_grid::UniformGrid2{T} # пространственная сетка
    time_grid::TimeGrid{T} # Временная сетка
    u::Array{Cdouble, 3} # массив решения
    f::Tuple{Matrix{Cdouble}, Matrix{Cdouble}} # массив потоков (`x` компоненты и `y` компоненты)


    function UniformGrid2Solution(spacial_grid::UniformGrid2{T}, time_grid::TimeGrid{T}) where T
        u = Array{Cdouble, 3}(undef, (spacial_grid.g[1].Ncells, spacial_grid.g[2].Ncells, time_grid.Ntime))
        f = (
            Matrix{Cdouble}(undef, (spacial_grid.g[1].Nfaces, spacial_grid.g[2].Ncells)),
            Matrix{Cdouble}(undef, (spacial_grid.g[1].Ncells, spacial_grid.g[2].Nfaces))
        )

        new{T}(spacial_grid, time_grid, u, f)
    end
end


"""
Задача Дирихле для двумерного квазилинейного уравнения теплопроводности
"""
struct HeatProblem2{K1<:Function, K2<:Function, F<:Function, Mu_1<:Function, Mu1<:Function, Mu_2<:Function, Mu2<:Function, U<:Function}
    k1::K1
    k2::K2
    f::F
    mu_1::Mu_1
    mu1::Mu1
    mu_2::Mu_2
    mu2::Mu2
    u0::U
end


"""
Обновление `f`-- потоков на временном слое `j`
"""
function update_flux!(
    problem::HeatProblem2, solution::UniformGrid2Solution,
    u_1::Vector{Cdouble}, u1::Vector{Cdouble},
    u_2::Vector{Cdouble}, u2::Vector{Cdouble}, j::Int
    )

    for k in 1:solution.spacial_grid.g[2].Ncells
        u_1[k] = 2 * problem.mu_1(
            solution.spacial_grid.g[2].cells[k],
            solution.time_grid.t[j]
        ) - solution.u[1, k, j]

        u1[k] = 2 * problem.mu1(
            solution.spacial_grid.g[2].cells[k],
            solution.time_grid.t[j]
        ) - solution.u[end, k, j]
    end

    for i in 1:solution.spacial_grid.g[1].Ncells
        u_2[i] = 2 * problem.mu_2(
            solution.spacial_grid.g[1].cells[i],
            solution.time_grid.t[j]
        ) - solution.u[i, 1, j]

        u2[i] = 2 * problem.mu2(
            solution.spacial_grid.g[1].cells[i],
            solution.time_grid.t[j]
        ) - solution.u[i, end, j]
    end

    for k in 1:solution.spacial_grid.g[2].Ncells
        solution.f[1][1, k] = -problem.k1(
            0.5 * (u_1[k] + solution.u[1, k, j])
        ) * (
            solution.u[1, k, j] - u_1[k]
        ) / solution.spacial_grid.g[1].h

        solution.f[1][end, k] = -problem.k1(
            0.5 * (solution.u[end, k, j] + u1[k])
        ) * (
            u1[k] - solution.u[end, k, j]
        ) / solution.spacial_grid.g[1].h
    end

    for i in 1:solution.spacial_grid.g[1].Ncells
        solution.f[2][i, 1] = -problem.k2(
            0.5 * (u_2[i] + solution.u[i, 1, j])
        ) * (
            solution.u[i, 1, j] - u_2[i]
        ) / solution.spacial_grid.g[2].h

        solution.f[2][i, end] = -problem.k2(
            0.5 * (solution.u[i, end, j] + u2[i])
        ) * (
            u2[i] - solution.u[i, end, j]
        ) / solution.spacial_grid.g[2].h
    end

    for i in 2:(solution.spacial_grid.g[1].Nfaces - 1)
        for k in 2:(solution.spacial_grid.g[2].Nfaces - 1)
            solution.f[1][i, k] = -problem.k1(
                0.5 * (solution.u[i - 1, k, j] + solution.u[i, k, j])
            ) * (
                solution.u[i, k, j] - solution.u[i - 1, k, j]
            ) / solution.spacial_grid.g[1].h

            solution.f[2][i, k] = -problem.k2(
                0.5 * (solution.u[i, k - 1, j] + solution.u[i, k, j])
            ) * (
                solution.u[i, k, j] - solution.u[i, k - 1, j]
            ) / solution.spacial_grid.g[2].h
        end
    end
end


"""
Учёт начальных условий
"""
function initial_conditions!(problem::HeatProblem2, solution::UniformGrid2Solution, u_1::Vector{Cdouble}, u1::Vector{Cdouble}, u_2::Vector{Cdouble}, u2::Vector{Cdouble})
    for i in 1:solution.spacial_grid.g[1].Ncells
        for k in 1:solution.spacial_grid.g[2].Ncells
            solution.u[i, k, 1] = problem.u0(
                solution.spacial_grid.g[1].cells[i],
                solution.spacial_grid.g[2].cells[k]
            )
        end
    end

    update_flux!(problem, solution, u_1, u1, u_2, u2, 1)
end


"""
Обновление `u` при временном шаге `j ↦ j + 1`
"""
function update_u!(
    problem::HeatProblem2, solution::UniformGrid2Solution, 
    u_1::Vector{Cdouble}, u1::Vector{Cdouble}, u_2::Vector{Cdouble}, u2::Vector{Cdouble},
    j::Int
    )

    for i = 1:(solution.spacial_grid.g[1].Ncells)
        for k = 1:(solution.spacial_grid.g[2].Ncells)
            solution.u[i, k, j + 1] = solution.u[i, k, j] - solution.time_grid.τ * (
                (solution.f[1][i + 1, k] - solution.f[1][i, k]) / solution.spacial_grid.g[1].h +
                (solution.f[2][i, k + 1] - solution.f[2][i, k]) / solution.spacial_grid.g[2].h
            ) + solution.time_grid.τ * problem.f(
                solution.spacial_grid.g[1].cells[i],
                solution.spacial_grid.g[2].cells[k],
                solution.time_grid.t[j]
            )
        end
    end
end


function update_u_loc!(
    problem::HeatProblem2, solution::UniformGrid2Solution, j::Int,
    α::Tuple{Vector{Cdouble}, Vector{Cdouble}},
    β::Tuple{Vector{Cdouble}, Vector{Cdouble}},
    temp_x::Tuple{Vector{Cdouble}, Vector{Cdouble}},
    A::Tuple{Vector{Cdouble}, Vector{Cdouble}},
    B::Tuple{Vector{Cdouble}, Vector{Cdouble}},
    C::Tuple{Vector{Cdouble}, Vector{Cdouble}},
    D::Tuple{Vector{Cdouble}, Vector{Cdouble}}
    )

    k1 = problem.k1
    k2 = problem.k2
    f = problem.f
    mu_1 = problem.mu_1
    mu1 = problem.mu1
    mu_2 = problem.mu_2
    mu2 = problem.mu2

    u = solution.u
    x = solution.spacial_grid.g[1].cells
    y = solution.spacial_grid.g[2].cells
    hx = solution.spacial_grid.g[1].h
    hy = solution.spacial_grid.g[2].h
    Nx = solution.spacial_grid.g[1].Ncells
    Ny = solution.spacial_grid.g[2].Ncells
    t = solution.time_grid.t
    τ = solution.time_grid.τ

    w = zeros(Nx, Ny)

    t_half = t[j] + 0.5 * τ

    # Сначала решаем задачу вдоль оси x для всех k
    for k = 1:Ny

        # Создаём трёхдиагональную матрицу
        for i = 2:(Nx - 1)
            A[1][i - 1] = -(1 / hx^2) * k1(0.5 * (u[i - 1, k, j] + u[i, k, j]))
            B[1][i] = 1 / τ + (1 / hx^2) * (
                k1(0.5 * (u[i - 1, k, j] + u[i, k, j])) +
                k1(0.5 * (u[i + 1, k, j] + u[i, k, j]))
            )
            C[1][i] = -(1 / hx^2) * k1(0.5 * (u[i + 1, k, j] + u[i, k, j]))
            D[1][i] = (1 / τ) * u[i, k, j] + 0.5 * f(x[i], y[k], t[j])
        end

        # учёт граничных условий
        u_1 = 2 * mu_1(y[k], t[j]) - u[1, k, j]
        u1 = 2 * mu1(y[k], t[j]) - u[end, k, j]

        B[1][1] = hx^2 + τ * k1(0.5 * (u[1, k, j] + u[2, k, j])) + 2 * τ * k1(0.5 * (u_1 + u[1, k, j]))
        C[1][1] = -τ * k1(0.5 * (u[1, k, j] + u[2, k, j]))
        D[1][1] = τ * hx^2 * f(x[1], y[k], t[j]) + hx^2 * u[1, k, j] + 2 * mu_1(y[k], t_half) * τ * k1(0.5 * (u_1 + u[1, k, j]))

        A[1][end] = -τ * k1(0.5 * (u[end - 1, k, j] + u[end, k, j]))
        B[1][end] = hx^2 + τ * k1(0.5 * (u[end - 1, k, j] + u[end, k, j])) + 2 * τ * k1(0.5 * (u1 + u[end, k, j]))
        D[1][end] = τ * hx^2 * f(x[end], y[k], t[j]) + hx^2 * u[end, k, j] + 2 * mu1(y[k], t_half) * τ * k1(0.5 * (u[end, k, j] + u1))

        # получаем решение
        tridiagonal_matrix_algorithm!(α[1], β[1], temp_x[1], A[1], B[1], C[1], D[1])
        w[:, k] .= temp_x[1]
    end

    # Теперь то же самое вдоль оси y для всех i:
    for i = 1:Nx

        # Создаём трёхдиагональную матрицу
        for k = 2:(Ny - 1)
            A[2][k - 1] = -(1 / hy^2) * k2(0.5 * (w[i, k - 1] + w[i, k]))
            B[2][k] = 1 / τ + (1 / hy^2) * (
                k2(0.5 * (w[i, k] + w[i, k - 1])) +
                k2(0.5 * (w[i, k + 1] + w[i, k]))
            )
            C[2][k] = -(1 / hy^2) * k2(0.5 * (w[i, k] + w[i, k + 1]))
            D[2][k] = (1 / τ) * w[i, k] + 0.5 * f(x[i], y[k], t_half)
        end

        # учёт граничных условий
        u_2 = 2 * mu_2(x[i], t_half) - w[i, 1]
        u2 = 2 * mu2(x[i], t_half) - w[i, end]

        B[2][1] = hy^2 + τ * k2(0.5 * (w[i, 1] + w[i, 2])) + 2 * τ * k2(0.5 * (u_2 + w[i, 1]))
        C[2][1] = -τ * k2(0.5 * (w[i, 1] + w[i, 2]))
        D[2][1] = τ * hy^2 * f(x[i], y[1], t_half) + hy^2 * w[i, 1] + 2 * mu_2(x[i], t[j + 1]) * τ * k2(0.5 * (u_2 + w[i, 1]))

        A[2][end] = -τ * k2(0.5 * (w[i, end - 1] + w[i, end]))
        B[2][end] = hy^2 + τ * k2(0.5 * (w[i, end - 1] + w[i, end])) + 2 * τ * k2(0.5 * (u2 + w[i, end]))
        D[2][end] = τ * hy^2 * f(x[i], y[end], t_half) + hy^2 * w[i, end] + 2 * mu2(x[i], t[j + 1]) * τ * k2(0.5 * (w[i, end] + u2))

        # получаем решение 
        tridiagonal_matrix_algorithm!(α[2], β[2], temp_x[2], A[2], B[2], C[2], D[2])
        u[i, :, j + 1] .= temp_x[2]
    end

end

"""
Цикл по всем временным слоям (основной цикл разностной схемы)
"""
function time_layers_loop!(
    problem::HeatProblem2, solution::UniformGrid2Solution,
    u_1::Vector{Cdouble}, u1::Vector{Cdouble}, u_2::Vector{Cdouble}, u2::Vector{Cdouble}
    )

    # Векторы для прогонки
    A::Tuple{Vector{Float64}, Vector{Float64}} = (
        Vector{Float64}(undef, solution.spacial_grid.g[1].Ncells - 1),
        Vector{Float64}(undef, solution.spacial_grid.g[2].Ncells - 1)
    )
    B::Tuple{Vector{Float64}, Vector{Float64}} = (
        Vector{Float64}(undef, solution.spacial_grid.g[1].Ncells),
        Vector{Float64}(undef, solution.spacial_grid.g[2].Ncells)
    )
    C::Tuple{Vector{Float64}, Vector{Float64}} = (
        Vector{Float64}(undef, solution.spacial_grid.g[1].Ncells - 1),
        Vector{Float64}(undef, solution.spacial_grid.g[2].Ncells - 1)
    )
    D::Tuple{Vector{Float64}, Vector{Float64}} = (
        Vector{Float64}(undef, solution.spacial_grid.g[1].Ncells),
        Vector{Float64}(undef, solution.spacial_grid.g[2].Ncells)
    )

    # Временные векторы для прогонки
    α::Tuple{Vector{Float64}, Vector{Float64}} = (
        Vector{Float64}(undef, solution.spacial_grid.g[1].Ncells),
        Vector{Float64}(undef, solution.spacial_grid.g[2].Ncells)
    )
    β::Tuple{Vector{Float64}, Vector{Float64}} = (
        Vector{Float64}(undef, solution.spacial_grid.g[1].Ncells),
        Vector{Float64}(undef, solution.spacial_grid.g[2].Ncells)
    )
    temp_x::Tuple{Vector{Float64}, Vector{Float64}} = (
        Vector{Float64}(undef, solution.spacial_grid.g[1].Ncells),
        Vector{Float64}(undef, solution.spacial_grid.g[2].Ncells)
    )

    for j = 1:(solution.time_grid.Ntime - 1)
        # update_u!(problem, solution, u_1, u1, u_2, u2, j)
        update_u_loc!(problem, solution, j, α, β, temp_x, A, B, C, D)
        update_flux!(problem, solution, u_1, u1, u_2, u2, j + 1)
    end
end


"""
Решение проблемы `HeatProblem2`
"""
function solve_PDE(problem::HeatProblem2, spacial_grid::UniformGrid2, time_grid::TimeGrid)

    # Создание структуры, хранящей решение:
    solution = UniformGrid2Solution(spacial_grid, time_grid)

    # Вспомогательные векторы
    u_1 = Vector{Cdouble}(undef, solution.spacial_grid.g[2].Ncells)
    u1 = Vector{Cdouble}(undef, solution.spacial_grid.g[2].Ncells)
    u_2 = Vector{Cdouble}(undef, solution.spacial_grid.g[1].Ncells)
    u2 = Vector{Cdouble}(undef, solution.spacial_grid.g[1].Ncells)

    # Учёт начальных условий
    initial_conditions!(problem, solution, u_1, u1, u_2, u2)

    # Цикл по временным слоям
    time_layers_loop!(problem, solution, u_1, u1, u_2, u2)

    # Возвращаем решение
    return solution
end