include("../solver.jl")
include("../plotter.jl")


function test()
    println("Preparing...")
    k(u) = 1.0
    f(x, t) = 0.0
    mu_(t) = 0.0
    mu(t) = 0.0
    u0(x) = x * sin(3π * x)
    problem = HeatProblem(k, f, mu_, mu, u0)

    spacial_grid = UniformGrid(100)
    time_grid = TimeGrid(0.05, 1001)
    println("done!")

    println("Solving...")
    sol = solve_PDE(problem, spacial_grid, time_grid)
    println("done!")

    println("Plotting grid...")
    plot_grid(spacial_grid, filename="grid_1.pdf")

    println("Plotting solution...")
    plot_gif(sol)
    println("done!")
end