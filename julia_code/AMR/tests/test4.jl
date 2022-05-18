include("../solver.jl")
include("../solver2.jl")
include("../plotter.jl")


function test()
    println("Preparing...")
    k(u) = 1.0
    f(x, y, t) = 0.0
    mu_1(y, t) = 1.0
    mu1(y, t) = 0.0
    mu_2(x, t) = 1.0 - 0.5 * x
    mu2(x, t) = 1.0 - 0.5 * x
    u0(x, y) = 0.5
    problem = HeatProblem2(k, k, f, mu_1, mu1, mu_2, mu2, u0)

    spacial_grid = UniformGrid2((2.0, 1.0), (30, 30))
    time_grid = TimeGrid(0.1, 1001)

    println("done!")

    println("Plotting grid...")
    plot_grid(spacial_grid, filename=joinpath("output", "grid_test2.pdf"))
    println("done!")

    println("Solving...")
    @time sol = solve_PDE(problem, spacial_grid, time_grid)
    println("done!")

    println("Plotting...")
    plot_gif(sol; filename=joinpath("output", "solution_test4.gif"))
    println("done!")
end