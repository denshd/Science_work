include("grids.jl")
include("problem.jl")
include("plotter.jl")
include("tridiagonal_matrix_algorithm.jl")


function initial_conditions! end
function advance_level! end
function update_u! end
function time_step! end
function set_boundary_values end
function interpolate_fine_to_coarse end


function initial_conditions!(problem::HeatProblem2, block::Block{T}) where T
    for i in 1:block.spacial_grid.g[1].N
        for j in 1:block.spacial_grid.g[2].N
            block.u_old[i, j] = problem.u0(block.spacial_grid.g[1].x[i], block.spacial_grid.g[2].x[j])
        end
    end
    block.u_new .= block.u_old
end


function initial_conditions!(problem::HeatProblem2, L::Level{T}) where T
    
    # По всем блокам на данном уровне
    for i = 1:L.M
        initial_conditions!(problem, L.blocks[i])
    end

    # Если есть подуровень, делаем то же самое
    if has_sublevel(L)
        initial_conditions!(problem, L.sublevel)
    end
end


function initial_conditions!(problem::HeatProblem2, solution::Solution{T}) where T
    first_level = solution.levels[1]

    initial_conditions!(problem, first_level)
end



function solve_PDE(problem::HeatProblem2, solution::Solution{T}) where T
    
    # Учёт самых начальных условий
    initial_conditions!(problem, solution)

    # Идём по всем временным слоям:
    for j in 1:(solution.time_grid.N - 1)
        # Вот он новый уровень, его и будем "advance'ировать"
        solution.levels[j + 1] = deepcopy(solution.levels[j])

        advance_level!(problem, solution.levels[j + 1])
    end
end


function advance_level!(problem::HeatProblem2, L::Level{T}) where T

    for i in 1:get_level_refinement_ratio(L)
        update_u!(problem, L, i) # Обновили решение
        if (has_sublevel(L)) # Если есть подуровень, "адвансируем и его"
            advance_level!(problem, L.sublevel)
        end
        L.t_curr += L.Δt # Обновляем время (решение, полученное выше в update_u как раз и будет решением в этот новый момент времени t_curr)
    end
    interpolate_fine_to_coarse(L)
end


function set_physical_boundary(problem::HeatProblem2, level::Level{T}) where T
    @. level.blocks[1].u_new[1, :] = problem.mu_1(level.blocks[1].spacial_grid.g[2].x, level.t_curr + level.Δt)
    @. level.blocks[1].u_new[end, :] = problem.mu1(level.blocks[1].spacial_grid.g[2].x, level.t_curr + level.Δt)
    @. level.blocks[1].u_new[:, 1] = problem.mu_2(level.blocks[1].spacial_grid.g[1].x, level.t_curr + level.Δt)
    @. level.blocks[1].u_new[:, end] = problem.mu2(level.blocks[1].spacial_grid.g[1].x, level.t_curr + level.Δt)
end


# TODO: добавить нормальный подсчёт `w_boundary`
# TODO: разобраться в ошибке
function set_boundary_values(problem::HeatProblem2, level::Level{T}, time_i::Int) where T

    if level.level_number == 1
        set_physical_boundary(problem, level)
        return nothing
    end

    suplevel = get_suplevel(level)
    for block in level.blocks # задаём граничные условия для каждого блока на этом уровне
        supblock = suplevel.blocks[block.supblock_index] # надблок `supblock` для данного блока `block`

        # координаты блока в надблоке в надблочных индексах
        i1 = block.supblock_position[1][1] 
        i2 = block.supblock_position[1][2]
        j1 = block.supblock_position[2][1]
        j2 = block.supblock_position[2][2]

        # Точная пространственная интерполяция с линейной временной
        is = 0
        for i = 1:2:block.spacial_grid.g[1].N # идём по совпадающим точкам по оси x (для "верхней палки" и "нижней палки")
            block.u_new[i, 1] = (1 - 0.5 * time_i) * supblock.u_old[i1 + is, j1] + 0.5 * time_i * supblock.u_new[i1 + is, j1]
            block.u_new[i, end] = (1 - 0.5 * time_i) * supblock.u_old[i1 + is, j2] + 0.5 * time_i * supblock.u_new[i1 + is, j2]
            is += 1
        end
        js = 0
        for j = 1:2:block.spacial_grid.g[2].N # идём по совпадающим точкам по оси y (для "левой палки" и "правой палки")
            block.u_new[1, j] = (1 - 0.5 * time_i) * supblock.u_old[i1, j1 + js] + 0.5 * time_i * supblock.u_new[i1, j1 + js]
            block.u_new[end, j] = (1 - 0.5 * time_i) * supblock.u_old[i2, j1 + js] + 0.5 * time_i * supblock.u_new[i2, j1 + js]
            js += 1
        end

        # Линейная пространственная интерполяция
        for i = 2:2:block.spacial_grid.g[1].N
            block.u_new[i, 1] = 0.5 * (block.u_new[i - 1, 1] + block.u_new[i + 1, 1])
            block.u_new[i, end] = 0.5 * (block.u_new[i - 1, end] + block.u_new[i + 1, end])
        end
        for j = 2:2:block.spacial_grid.g[2].N
            block.u_new[1, j] = 0.5 * (block.u_new[1, j - 1] + block.u_new[1, j + 1])
            block.u_new[end, j] = 0.5 * (block.u_new[end, j - 1] + block.u_new[end, j + 1])
        end
    end
end


function interpolate_fine_to_coarse(level::Level{T}) where T
    if level.level_number == 1
        return nothing
    end
    suplevel = get_suplevel(level)
    for block in level.blocks # интерполируем для каждого блока на этом уровне
        supblock = suplevel.blocks[block.supblock_index] # надблок `supblock` для данного блока `block`

        # координаты блока в надблоке в надблочных индексах
        i1 = block.supblock_position[1][1] 
        i2 = block.supblock_position[1][2]
        j1 = block.supblock_position[2][1]
        j2 = block.supblock_position[2][1]

        # Точная пространственная интерполяция с линейной временной
        is = 0
        js = 0
        for i = i1:i2
            for j = j1:j2
                supblock.u_new[i, j] = supblock.u_new[i, j]
                supblock.u_new[i, j] = block.u_new[2*(i - i1 + 1) - 1, 2*(j - j1 + 1) - 1]
            end
        end
    end

end

function update_u!(problem::HeatProblem2, L::Level{T}, time_i::Int) where T

    # Обновляем то, что было новое, теперь старое, то, что было старое, уже и не нужно
    for i in 1:L.M
        L.blocks[i].u_old .= L.blocks[i].u_new
    end

    # Установка граничных условий (тут учитываются и физические, и coarse-fine, и так далее)
    set_boundary_values(problem, L, time_i)

    # Идём по всем блокам в уровне L
    for i in 1:L.M
        time_step!(problem, L.blocks[i], L.t_curr, L.Δt)
        
        # Вот здесь будет обновление потоков, если это действительно будет необходимо (вряд ли, ибо мы считаем неявной схемой)
    end
end


function time_step!(problem::HeatProblem2, block::Block{T}, t_curr::Cdouble, Δt::Cdouble) where T
    A = block.memory_for_time_step.A
    B = block.memory_for_time_step.B
    C = block.memory_for_time_step.C
    D = block.memory_for_time_step.D
    α = block.memory_for_time_step.α
    β = block.memory_for_time_step.β
    temp_x = block.memory_for_time_step.temp_x

    u_new = block.u_new
    u_old = block.u_old
    w_boundary = block.w_boundary
    w = block.w

    hx = block.spacial_grid.g[1].h
    hy = block.spacial_grid.g[2].h
    Nx = block.spacial_grid.g[1].N
    Ny = block.spacial_grid.g[2].N
    x = block.spacial_grid.g[1].x
    y = block.spacial_grid.g[2].x

    k1 = problem.k1
    k2 = problem.k2
    mu_1 = problem.mu_1
    mu1 = problem.mu1
    mu_2 = problem.mu_2
    mu2 = problem.mu2
    f = problem.f

    t_half = t_curr + 0.5 * Δt

    # Здесь Локально-одномерная схема

    # Сначала решаем задачу вдоль оси x (для всех y[k])
    for k = 1:Ny

        # Создаём трёхдиагональную матрицу
        for i = 2:(Nx - 1)
            A[1][i - 1] = -(1 / hx^2) * k1(0.5 * (u_old[i - 1, k] + u_old[i, k]))
            B[1][i] = 1 / Δt + (1 / hx^2) * (
                k1(0.5 * (u_old[i - 1, k] + u_old[i, k])) +
                k1(0.5 * (u_old[i + 1, k] + u_old[i, k]))
            )
            C[1][i] = -(1 / hx^2) * k1(0.5 * (u_old[i + 1, k] + u_old[i, k]))
            D[1][i] = (1 / Δt) * u_old[i, k] + 0.5 * f(x[i], y[k], t_curr)
        end
        # Граничные значения уже были учтены в set_boundary_values()
        B[1][1] = B[1][end] = 1
        C[1][1] = A[1][end] = 0
        D[1][1] = u_new[1, k]
        D[1][end] = u_new[end, k]

        # Получаем временный вектор w
        tridiagonal_matrix_algorithm!(α[1], β[1], temp_x[1], A[1], B[1], C[1], D[1])
        w[:, k] .= temp_x[1]
    end

    # Затем для вдоль оси y (для всех x[i])
    for i = 1:Nx

        # Создаём трёхдиагональную матрицу
        for k = 2:(Ny - 1)
            A[2][k - 1] = -(1 / hy^2) * k2(0.5 * (w[i, k] + w[i, k - 1]))
            B[2][k] = 1 / Δt + (1 / hy^2) * (
                k2(0.5 * (w[i, k] + w[i, k - 1])) +
                k2(0.5 * (w[i, k + 1] + w[i, k]))
            )
            C[2][k] = -(1 / hy^2) * k2(0.5 * (w[i, k] + w[i, k + 1]))
            D[2][k] = (1 / Δt) * w[i, k] + 0.5 * f(x[i], y[k], t_half)
        end
        # Опять же, граничные значения уже были учтены в set_boundary_values()
        B[2][1] = B[2][end] = 1
        C[2][1] = A[2][end] = 0
        D[2][1] = u_new[i, 1]
        D[2][end] = u_new[i, end]

        # Получаем решение
        tridiagonal_matrix_algorithm!(α[2], β[2], temp_x[2], A[2], B[2], C[2], D[2])
        u_new[i, :] .= temp_x[2]
    end
end