function example_2()
    f(x, t, u) = 0.0
    mu1(t) = 10 * √t
    mu2(t) = 0.0
    ϕ(x) = x < 0.5 ? (2 * √(5 * (0.5 - x))) : 0.0
    k(u) = 0.5 * u^2

    sys = Heat_quazilinear_1(k, f, mu1, mu2, ϕ)

    # Создаём интересную сетку
    x1 = 0
    x2 = 1
    t1 = 0.1
    t2 = 0.2
    Nx = 11
    Nt = 101
    g = GeneralGrid(
        sort(unique(vcat(
            collect(range(start=0, stop=1, length=Nx)),
            collect(range(start=0.5, stop=1, length=2 * Nx - 1)),
            collect(range(start=0.75, stop=1, length=2 * Nx - 1)),
            collect(range(start=0.825, stop=1, length=2 * Nx - 1))
        ))),
        collect(range(t1, t2, length=Nt))
    )

    # Решаем
    sol = solve_PDE(sys, g)
    @time sol = solve_PDE(sys, g)

    plot_grid(g, filename="output/grid.pdf")
    plot_gif(sol, filename="output/example_2.gif")
end