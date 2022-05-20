include("../solver.jl")
include("../solver2.jl")
include("../plotter.jl")


function test()
    println("Preparing...")
    k1(u) = u^4
    k2(u) = u^4
    f(x, y, t) = 3 * ℯ^(-((x - 0.25)^2 + (y - 0.25)^2) / 0.1^2) * abs(sin(100 * t))
    mu(x, t) = 0.0

    u0(x, y) = 3 * ℯ^(-((x - 0.5)^2 + (y - 0.5)^2) / 0.1^2)
    problem = HeatProblem2(k1, k2, f, mu, mu, mu, mu, u0)

    L1 = 1
    L2 = 1
    N = (200, 200)

    T = 0.5
    Nt = 1000

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
    plot_gif(sol; filename=joinpath("output", "solution_test6.gif"))
    println("done!")
end