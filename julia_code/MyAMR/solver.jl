include("grids.jl")
include("problem.jl")
include("plotter.jl")


function initial_conditions! end
function advance_level! end
function update_u! end
function time_step! end



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

    # В зависимости от того, насколько "глубокий" этот уровень, столько раз и обновляем его
    for i in 1:(2^(L.level_number - 1))
        # Вот здесь я думаю, перед тем, как обновлять решение, нужно интерполировать все граничные/начальные условия... TODO: разобраться
        update_u!(problem, L) # Обновили решение
        if (has_sublevel(L)) # Если есть подуровень, "адвансируем и его"
            advance_level!(problem, L.sublevel)
        end
        L.t_curr += L.Δt # Обновляем время (решение, полученное выше в update_u как раз и будет решением в этот новый момент времени t_curr)
    end
end


function update_u!(problem::HeatProblem2, L::Level{T}) where T

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

    hx = block.spacial_grid.g[1].h
    hy = block.spacial_grid.g[2].h
    Nx = block.spacial_grid.g[1].N
    Ny = block.spacial_grid.g[2].N

    k1 = problem.k1
    k2 = problem.k2
    mu_1 = problem.mu_1
    mu1 = problem.mu1
    mu_2 = problem.mu_2
    mu2 = problem.mu2
    f = problem.f

    t_half = t_curr + 0.5 * Δt

    # Интерполяция 
    # Здесь Локально-одномерная схема

    
end