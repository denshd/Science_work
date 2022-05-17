"""
Равномерная сетка из клеток и граней
"""
struct UniformGrid{T <: AbstractRange{Cdouble}}
    Ncells::Int # число клеточек
    Nfaces::Int # число палочек
    h::Cdouble # шаг сетки
    half_h::Cdouble # h / 2
    cells::T # координаты клеточек в виде Range
    faces::T # координаты палочек в виде Range


    """
    Создаёт равномерную сетку на `[0, 1]` с числом клеток `Ncells`
    """
    function UniformGrid(Ncells::Integer)
        Nfaces = Ncells + 1
        faces = range(start=0.0, stop=1, length=Nfaces)
        cells = range(start=(0.0 + step(faces)), step=step(faces), length=Ncells)

        new{typeof(cells)}(Ncells, Nfaces, step(cells), 0.5 * step(cells), cells, faces)
    end


    """
    Создаёт равномерную сетку на `[0, 1]` с шагом `h`
    """
    function UniformGrid(h::Real)
        faces = range(start=0, stop=1, step=h)
        cells = range(start=(0 + h), step=h, length=(length(faces) - 1))

        new{typeof(cells)}(length(cells), length(faces), h, 0.5 * h, cells, faces)
    end


    """
    Создаёт равномерную сетку на `[0, L]` с числом клеток `Ncells`
    """
    function UniformGrid(L::Real, Ncels::Integer)
        faces = range(start=0, stop=L, length=(Ncells + 1))
        cells = range(start=(0 + step(faces)), step=step(faces), length=Ncells)

        new{typeof(cells)}(length(cells), length(faces), step(cells), 0.5 * step(cells), cells, faces)
    end


    """
    Создаёт равномерную сетку на `[0, L]` с шагом `h`
    """
    function UniformGrid(L::Real, h::Real)
        faces = range(start=0, stop=L, step=h)
        cells = range(start=0, step=h, length=(length(faces) - 1))

        new{typeof(cells)}(length(cells), length(faces), h, 0.5 * h, cells, faces)
    end
end


"""
Равномерная временная сетка
"""
struct TimeGrid{T <: AbstractRange{Cdouble}}
    Ntime::Int # число временных точек
    τ::Cdouble # шаг времени
    τ_half::Cdouble # τ / 2
    t::T # сами времена в виде Range


    """
    Создаёт равномерную временную сетку на `[0, T]` с числом точек Ntime
    """
    function TimeGrid(T::Real, Ntime::Integer)
        t = range(start=0, stop=T, length=Ntime)

        new{typeof(t)}(Ntime, step(t), 0.5 * step(t), t)
    end


    """
    Создаёт равномерную временную сетку на `[t1, t2]` с числом точек Ntime
    """
    function TimeGrid((t1, t2)::Tuple{Real, Real}, Ntime::Integer)
        t = range(start=t1, stop=t2, length=Ntime)
        new{typeof(t)}(Ntime, step(t), 0.5 * step(t), t)
    end
end


"""
Представление решения
"""
struct UniformGridSolution{T <: AbstractRange{Cdouble}}
    spacial_grid::UniformGrid{T} # пространственная сетка
    time_grid::TimeGrid{T} # временная сетка
    u::Matrix{Cdouble} # матрица решения
    f::Vector{Cdouble} # вектор потоков


    function UniformGridSolution(spacial_grid::UniformGrid{T}, time_grid::TimeGrid{T}) where T
        u = Matrix{Cdouble}(undef, (spacial_grid.Ncells, time_grid.Ntime))
        f = Vector{Cdouble}(undef, spacial_grid.Nfaces)
        new{T}(spacial_grid, time_grid, u, f)
    end
end


"""
Задача Дирихле для одномерного квазилинейного уравнения теплопроводности
"""
struct HeatProblem{K<:Function, F<:Function, Mu_<:Function, Mu<:Function, U<:Function}
    k::K
    f::F
    mu_::Mu_
    mu::Mu
    u0::U
end


"""
Учёт начальных условий
"""
function initial_conditions(problem::HeatProblem, solution::UniformGridSolution)
    @. solution.u[:, 1] = problem.u0(solution.spacial_grid.cells)

    u_ = 2 * problem.mu_(solution.time_grid.t[1]) - solution.u[1, 1]
    u = 2 * problem.mu(solution.time_grid.t[1]) - solution.u[end, 1]
    solution.f[1] = -problem.k(0.5 * (u_ + solution.u[1, 1])) * (solution.u[1, 1] - u_) / solution.spacial_grid.h
    solution.f[end] = -problem.k(0.5 * (solution.u[end, 1] + u)) * (u - solution.u[end, 1]) / solution.spacial_grid.h
    for i = 2:(solution.spacial_grid.Nfaces - 1)
        solution.f[i] =  -problem.k(0.5 * (solution.u[i - 1, 1] + solution.u[i, 1])) * (solution.u[i, 1] - solution.u[i - 1, 1]) / solution.spacial_grid.h
    end
end


"""
Обновление `u` при временном шаге `j ↦ j + 1`
"""
function update_u(problem::HeatProblem, solution::UniformGridSolution, j::Int)
    for i = 1:(solution.spacial_grid.Ncells)
        solution.u[i, j + 1] = solution.u[i, j] - solution.time_grid.τ / solution.spacial_grid.h * (
            solution.f[i + 1] - solution.f[i]
        ) + solution.time_grid.τ * problem.f(solution.spacial_grid.cells[i], solution.time_grid.t[j])
    end
end


"""
Обновление `f`-- потоков на временном слое `j`
"""
function update_flux(problem::HeatProblem, solution::UniformGridSolution, j::Int)
    u_ = 2 * problem.mu_(solution.time_grid.t[j]) - solution.u[1, j]
    u = 2 * problem.mu(solution.time_grid.t[j]) - solution.u[end, j]
    solution.f[1] = -problem.k(0.5 * (u_ + solution.u[1, j])) * (solution.u[1, j] - u_) / solution.spacial_grid.h
    solution.f[end] = -problem.k(0.5 * (solution.u[end, j] + u)) * (u - solution.u[end, j]) / solution.spacial_grid.h
    for i = 2:(solution.spacial_grid.Nfaces - 1)
        solution.f[i] =  -problem.k(0.5 * (solution.u[i - 1, j] + solution.u[i, j])) * (solution.u[i, j] - solution.u[i - 1, j]) / solution.spacial_grid.h
    end
end


"""
Цикл по всем временным слоям (основной цикл разностной схемы)
"""
function time_layers_loop(problem::HeatProblem, solution::UniformGridSolution)
    for j = 1:(solution.time_grid.Ntime - 1)
        update_u(problem, solution, j)
        update_flux(problem, solution, j + 1)
    end
end


"""
Решение проблемы `HeatProblem`
"""
function solve_PDE(problem::HeatProblem, spacial_grid::UniformGrid, time_grid::TimeGrid)
    
    # Создание структуры, хранящей решение:
    solution = UniformGridSolution(spacial_grid, time_grid)

    # Учёт начальных условий
    initial_conditions(problem, solution)

    # Цикл по временным слоям
    time_layers_loop(problem, solution)

    # Возвращаем решение
    return solution
end