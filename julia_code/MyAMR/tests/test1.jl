include("../grids.jl")
include("../problem.jl")
include("../plotter.jl")
include("../solver.jl")


function test()
    println("Preparing...")
    k1(u) = u^4
    k2(u) = u^4
    f(x, y, t) = 3 * ℯ^(-((x - 0.25)^2 + (y - 0.25)^2) / 0.1^2) * abs(sin(100 * t))
    mu(x, t) = 0.0

    u0(x, y) = 3 * ℯ^(-((x - 0.5)^2 + (y - 0.5)^2) / 0.1^2)
    problem = HeatProblem2(k1, k2, f, mu, mu, mu, mu, u0)
    println("done!")

    # Здесь долгое создавание сетки...
    println("Creating grid...")
    solution = create_grid()
    println("done!")

    plot_grid(solution.levels[1], filename=joinpath("output", "first_result.pdf"))

    println("Solving...")
    solve_PDE(problem, solution)
    println("done!")

    println("Plotting solution...")
    plot_gif(solution, filename=joinpath("output", "first_result_gifka.gif"))
    println("done!")
    return nothing
end


function create_grid()

    L1 = 1
    L2 = 1
    N = (100, 100)

    T = 0.1
    Nt = 101
    time_grid = UniformGrid(T, Nt)
    Δt = time_grid.h

    main_block = create_block(
        UniformGrid2((101, 101)),
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

    # display(first_level)
    # if has_sublevel(first_level)
    #     display("What the hell?")
    # end
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