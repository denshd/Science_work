using Solver

function example_1()
    f(x, t, u) = 0.0
    mu1(t) = 10 * √t
    mu2(t) = 0.0
    u0(x) = x < 0.5 ? (2 * √(5 * (0.5 - x))) : 0.0
    k(u) = 0.5 * u^2

    sys = Heat_quazilinear_1(k, f, mu1, mu2, u0)

    x1 = 0
    x2 = 1
    t1 = 0.1
    t2 = 0.2
    Nx = 501
    Nt = 5001
    g = UniformGrid((x1, x2), (t1, t2), Nx, Nt)

    sol = solve_PDE(sys, g)
    @time sol = solve_PDE(sys, g)

    plot_gif(sol, filename="output/example_1.gif")
end