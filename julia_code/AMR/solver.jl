"""
Метод прогонки для решения трёхдиагональных систем линейных уравнений.
"""
function tridiagonal_matrix_algorithm(a::Vector{<:Number}, B::Vector{<:Number}, c::Vector{<:Number}, D::Vector{<:Number})::Vector{Float64}

    # размер матрицы
    n = length(B)

    # удобная проверка
    if (length(a) != length(B))
        A = vcat([0], copy(a))
    else
        A = copy(a)
    end
    if (length(c) != length(B))
        C = vcat(copy(c), [0])
    else
        C = copy(c)
    end

    # вспомогательный вектор y и коэффициенты алгоритма прогонки α и β, и столбец решение x
    # y = Vector{Float64}(undef, n)
    α = Vector{Float64}(undef, n)
    β = Vector{Float64}(undef, n)
    x = Vector{Float64}(undef, n)

    # Алгоритм прогонки:

    # прямой ход
    α[1] = - C[1] / B[1]
    β[1] = D[1] / B[1]
    for i=2:n-1
        α[i] = - C[i] / (B[i] + A[i] * α[i-1])
        β[i] = (D[i] - A[i] * β[i-1]) / (B[i] + A[i] * α[i-1])
    end
    β[n] = (D[n] - A[n] * β[n-1]) / (B[n] + A[n] * α[n-1])

    # обратный ход
    x[n] = β[n]
    for i = n-1:-1:1
        x[i] = α[i]*x[i+1] + β[i]
    end

    return x
end;


"""
Оптимальный по памяти алгоритм прогонки в случае, если его нужно запускать много раз
"""
function tridiagonal_matrix_algorithm!(
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
    end;
end


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
        cells = range(start=(0.0 + 0.5 * step(faces)), step=step(faces), length=Ncells)

        new{typeof(cells)}(Ncells, Nfaces, step(cells), 0.5 * step(cells), cells, faces)
    end


    """
    Создаёт равномерную сетку на `[0, 1]` с шагом `h`
    """
    function UniformGrid(h::Real)
        faces = range(start=0, stop=1, step=h)
        cells = range(start=(0 + 0.5 * h), step=h, length=(length(faces) - 1))

        new{typeof(cells)}(length(cells), length(faces), h, 0.5 * h, cells, faces)
    end


    """
    Создаёт равномерную сетку на `[0, L]` с числом клеток `Ncells`
    """
    function UniformGrid(L::Real, Ncells::Integer)
        faces = range(start=0, stop=L, length=(Ncells + 1))
        cells = range(start=(0 + 0.5 * step(faces)), step=step(faces), length=Ncells)

        new{typeof(cells)}(length(cells), length(faces), step(cells), 0.5 * step(cells), cells, faces)
    end


    """
    Создаёт равномерную сетку на `[0, L]` с шагом `h`
    """
    function UniformGrid(L::Real, h::Real)
        faces = range(start=0, stop=L, step=h)
        cells = range(start=(0 + 0.5 * step(faces)), step=h, length=(length(faces) - 1))

        new{typeof(cells)}(length(cells), length(faces), h, 0.5 * h, cells, faces)
    end


    """
    Создаёт равномерную сетку на `[x_b[1], x_b[2]]` с числом клеток `Ncells`
    """
    function UniformGrid(x_b::Tuple{Real, Real}, Ncells::Integer)
        faces = range(start=x_b[1], stop=x_b[2], length=(Ncells + 1))
        cells = range(start=(x_b[1] + 0.5 * step(faces)), step=step(faces), length=Ncells)

        new{typeof(cells)}(length(cells), length(faces), step(cells), 0.5 * step(cells), cells, faces)
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
Обновление `f`-- потоков на временном слое `j`
"""
function update_flux!(problem::HeatProblem, solution::UniformGridSolution, j::Int)
    u_ = 2 * problem.mu_(solution.time_grid.t[j]) - solution.u[1, j]
    u = 2 * problem.mu(solution.time_grid.t[j]) - solution.u[end, j]
    solution.f[1] = -problem.k(0.5 * (u_ + solution.u[1, j])) * (solution.u[1, j] - u_) / solution.spacial_grid.h
    solution.f[end] = -problem.k(0.5 * (solution.u[end, j] + u)) * (u - solution.u[end, j]) / solution.spacial_grid.h
    for i = 2:(solution.spacial_grid.Nfaces - 1)
        solution.f[i] =  -problem.k(0.5 * (solution.u[i - 1, j] + solution.u[i, j])) * (solution.u[i, j] - solution.u[i - 1, j]) / solution.spacial_grid.h
    end
end


"""
Учёт начальных условий
"""
function initial_conditions!(problem::HeatProblem, solution::UniformGridSolution)
    @. solution.u[:, 1] = problem.u0(solution.spacial_grid.cells)

    update_flux!(problem, solution, 1)
end


"""
Обновление `u` при временном шаге `j ↦ j + 1`
"""
function update_u!(problem::HeatProblem, solution::UniformGridSolution, j::Int)
    for i = 1:(solution.spacial_grid.Ncells)
        solution.u[i, j + 1] = solution.u[i, j] - solution.time_grid.τ / solution.spacial_grid.h * (
            solution.f[i + 1] - solution.f[i]
        ) + solution.time_grid.τ * problem.f(solution.spacial_grid.cells[i], solution.time_grid.t[j])
    end
end


"""
Обновление `u` при временном шаге `j ↦ j + 1` с помощью неявной схемы
"""
function update_u_loc!(
    problem::HeatProblem, solution::UniformGridSolution, j::Int,
    α::Vector{Cdouble}, β::Vector{Cdouble}, temp_x::Vector{Cdouble},
    A::Vector{Cdouble}, B::Vector{Cdouble}, C::Vector{Cdouble}, D::Vector{Cdouble}
    )

    # установка коэффициентов
    for i = 2:(solution.spacial_grid.Ncells - 1)
        A[i - 1] = -solution.time_grid.τ * problem.k(0.5 * (solution.u[i - 1, j] + solution.u[i, j]))
        B[i] = solution.spacial_grid.h^2 + solution.time_grid.τ * (
            problem.k(0.5 * (solution.u[i, j] + solution.u[i + 1, j])) +
            problem.k(0.5 * (solution.u[i - 1, j] + solution.u[i, j]))
        )
        C[i] = -solution.time_grid.τ * problem.k(0.5 * (solution.u[i, j] + solution.u[i + 1, j]))
        D[i] = solution.time_grid.τ * solution.spacial_grid.h^2 * problem.f(solution.spacial_grid.cells[i], solution.time_grid.t[j]) + solution.spacial_grid.h^2 * solution.u[i, j]
    end

    # учёт граничных условий
    u_ = 2 * problem.mu_(solution.time_grid.t[j]) - solution.u[1, j]
    u = 2 * problem.mu(solution.time_grid.t[j]) - solution.u[end, j]

    B[1] = solution.spacial_grid.h^2 + solution.time_grid.τ * problem.k(0.5 * (solution.u[1, j] + solution.u[2, j])) + 2 * solution.time_grid.τ * problem.k(0.5 * (u_ + solution.u[1, j]))
    C[1] = -solution.time_grid.τ * problem.k(0.5 * (solution.u[1, j] + solution.u[2, j]))
    D[1] = solution.time_grid.τ * solution.spacial_grid.h^2 * problem.f(solution.spacial_grid.cells[1], solution.time_grid.t[j]) + solution.spacial_grid.h^2 * solution.u[1, j] + 2 * problem.mu_(solution.time_grid.t[j + 1]) * solution.time_grid.τ * problem.k(0.5 * (u_ + solution.u[1, j]))

    A[end] = -solution.time_grid.τ * problem.k(0.5 * (solution.u[end - 1, j] + solution.u[end, j]))
    B[end] = solution.spacial_grid.h^2 + solution.time_grid.τ * problem.k(0.5 * (solution.u[end - 1, j] + solution.u[end, j])) + 2 * solution.time_grid.τ * problem.k(0.5 * (u + solution.u[end, j]))
    D[end] = solution.time_grid.τ * solution.spacial_grid.h^2 * problem.f(solution.spacial_grid.cells[end], solution.time_grid.t[j]) + solution.spacial_grid.h^2 * solution.u[end, j] + 2 * problem.mu(solution.time_grid.t[j + 1]) * solution.time_grid.τ * problem.k(0.5 * (solution.u[end, j] + u))

    tridiagonal_matrix_algorithm!(α, β, temp_x, A, B, C, D)
    solution.u[:, j + 1] .= temp_x
end


"""
Цикл по всем временным слоям (основной цикл разностной схемы)
"""
function time_layers_loop!(problem::HeatProblem, solution::UniformGridSolution)

    # Векторы для прогонки
    A = zeros(solution.spacial_grid.Ncells - 1)
    B = zeros(solution.spacial_grid.Ncells)
    C = zeros(solution.spacial_grid.Ncells - 1)
    D = zeros(solution.spacial_grid.Ncells)

    # Вспомогательные векторы для прогонки
    α = zeros(solution.spacial_grid.Ncells)
    β = zeros(solution.spacial_grid.Ncells)
    temp_x = zeros(solution.spacial_grid.Ncells)

    for j = 1:(solution.time_grid.Ntime - 1)
        # update_u!(problem, solution, j)
        update_u_loc!(problem, solution, j, α, β, temp_x, A, B, C, D)
        update_flux!(problem, solution, j + 1)
    end
end


"""
Решение проблемы `HeatProblem`
"""
function solve_PDE(problem::HeatProblem, spacial_grid::UniformGrid, time_grid::TimeGrid)
    
    # Создание структуры, хранящей решение:
    solution = UniformGridSolution(spacial_grid, time_grid)

    # Учёт начальных условий
    initial_conditions!(problem, solution)

    # Цикл по временным слоям
    time_layers_loop!(problem, solution)

    # Возвращаем решение
    return solution
end