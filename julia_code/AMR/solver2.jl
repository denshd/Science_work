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
    # @time begin
    # @. solution.f[1][1, :] = -problem.k1(
    #     0.5 * (u_1 + solution.u[1, :, j])
    # ) * (
    #     solution.u[1, :, j] - u_1
    # ) / solution.spacial_grid.g[1].h

    # @. solution.f[1][end, :] = -problem.k1(
    #     0.5 * (solution.u[end, :, j] + u1)
    # ) * (
    #     u1 - solution.u[end, :, j]
    # ) / solution.spacial_grid.g[1].h

    # @. solution.f[2][:, 1] = -problem.k2(
    #     0.5 * (u_2 + solution.u[:, 1, j])
    # ) * (
    #     solution.u[:, 1, j] - u_2
    # ) / solution.spacial_grid.g[2].h

    # @. solution.f[2][:, end] = -problem.k2(
    #     0.5 * (solution.u[:, end, j] + u2)
    # ) * (
    #     u2 - solution.u[:, end, j]
    # ) / solution.spacial_grid.g[2].h
    # end

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


"""
Цикл по всем временным слоям (основной цикл разностной схемы)
"""
function time_layers_loop!(
    problem::HeatProblem2, solution::UniformGrid2Solution,
    u_1::Vector{Cdouble}, u1::Vector{Cdouble}, u_2::Vector{Cdouble}, u2::Vector{Cdouble}
    )

    for j = 1:(solution.time_grid.Ntime - 1)
        update_u!(problem, solution, u_1, u1, u_2, u2, j)
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