include("../grids.jl")
include("../problem.jl")
include("../plotter.jl")
include("../solver.jl")


function test()
    println("Preparing...")
    k1(u) = 4 * u^4
    k2(u) = 0.25 * u^2
    f(x, y, t) = 0.0

    function mu_1(x2, t)
        if (t > (2*x2))
            return 0.5 * √(-1 + √(1 + 16 * (t - 2*x2)))
        else
            return 0
        end
    end

    function mu1(x2, t)
        if (t > (30 + 2*x2))
            return 0.5 * √(-1 + √(1 + 16 * (t - 30 - 2*x2)))
        else
            return 0
        end
    end

    function mu_2(x1, t)
        if (t > x1)
            return 0.5 * √(-1 + √(1 + 16 * (t - x1)))
        else
            return 0
        end
    end

    function mu2(x1, t)
        if (t > (x1 + 40))
            return 0.5 * √(-1 + √(1 + 16 * (t - x1 - 40)))
        else
            return 0
        end
    end
    u0(x, y) = 0.0
    problem = HeatProblem2(k1, k2, f, mu_1, mu1, mu_2, mu2, u0)

    println("Creating grid...")
    solution = create_grid()
    println("done!")

    plot_grid(solution.levels[1], filename=joinpath("output", "test2_grid.pdf"))

    println("Solving...")
    @time solve_PDE(problem, solution)
    println("done!")

    println("Plotting solution...")
    plot_gif(solution, filename=joinpath("output", "test2_solution.gif"))
    println("done!")
    return solution
end


function create_grid()

    L1 = 30
    L2 = 20
    N = (101, 101)

    T1 = 0
    T2 = 50
    Nt = 501
    time_grid = UniformGrid((T1, T2), Nt)
    Δt = time_grid.h

    main_block = create_block(
        UniformGrid2((L1, L2), N),
        1,
        1,
        ((1, 1), (1, 1))
    )

    second_block = create_subblock(
        main_block,
        1,
        ((21, 81), (21, 81))
    )

    third_blocks = [
        create_subblock(
            second_block,
            1,
            ((21, 41), (21, 41))
        )
        create_subblock(
            second_block,
            2,
            ((41, 81), (21, 41))
        )
        create_subblock(
            second_block,
            3,
            ((81, 101), (21, 41))
        )
        create_subblock(
            second_block,
            4,
            ((81, 101), (41, 81))
        )
        create_subblock(
            second_block,
            5,
            ((81, 101), (81, 111))
        )
        create_subblock(
            second_block,
            6,
            ((41, 81), (81, 111))
        )
        create_subblock(
            second_block,
            7,
            ((21, 41), (81, 111))
        )
        create_subblock(
            second_block,
            8,
            ((21, 41), (41, 81))
        )
    ]

    first_level = initialize_level(
        [main_block], 
        time_grid.x[1],
        Δt
    )

    add_sublevel!(first_level, [second_block])
    add_sublevel!(get_sublevel(first_level), third_blocks)

    levels = Vector{Level{get_type(first_level)}}(undef, Nt)
    levels[1] = first_level

    result = Solution(
        levels,
        time_grid
    )
    
    return result
end