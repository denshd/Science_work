include("../solver.jl")
include("../plotter.jl")


function test()
    println("Preparing...")
    k(u) = 0.5 * u^2
    f(x, t) = 0.0
    mu_(t) = 1 / √(0.9 - 8t)
    mu(t) = 0
    u0(x) = (x ≤ 0.5) ? (0.5 - x) / 0.15811388300841894 : 0
    problem = HeatProblem(k, f, mu_, mu, u0)

    spacial_grid = UniformGrid(100)
    time_grid = TimeGrid((0.1, 0.11), 50001)
    println("done!")

    println("Solving...")
    @time sol = solve_PDE(problem, spacial_grid, time_grid)
    println("done!")

    println("Plotting...")
    plot_gif(sol)
    println("done!")
end