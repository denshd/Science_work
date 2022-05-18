include("../solver.jl")
include("../solver2.jl")
include("../plotter.jl")


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

    L1 = 30
    L2 = 20
    N = (300, 200)

    T = 50
    Nt = 500

    spacial_grid = UniformGrid2((L1, L2), N)
    time_grid = TimeGrid(T, Nt)

    println("done!")

    # println("Plotting grid...")
    # plot_grid(spacial_grid, filename=joinpath("output", "2d_grid.pdf"))
    # println("done!")

    println("Solving...")
    @time sol = solve_PDE(problem, spacial_grid, time_grid)
    println("done!")

    println("Plotting...")
    plot_gif(sol; filename=joinpath("output", "solution_test5.gif"))
    println("done!")
end