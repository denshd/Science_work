using Solver

function example_3()
    # Неоднородность уравнения
    f(x, y, t) = 0.0

    # Граничные условия
    mu_1(y, t) = 1.0
    mu1(y, t) = 0.0
    mu_2(x, t) = 1.0 - (x / 2)
    mu2(x, t) = 1.0 - (x / 2)

    # Начальное услвоие
    u0(x, y) = 0.5

    # Задача
    sys = Heat2_linear_const_1(f, mu_1, mu1, mu_2, mu2, u0)

    x1 = 0
    x2 = 2

    y1 = 0
    y2 = 1

    t1 = 0
    t2 = 0.2

    Nx = 51
    Ny = 51
    Nt = 501

    # Сетка
    g = UniformGrid_2(((x1, x2), (y1, y2)), (t1, t2), (Nx, Ny), Nt)

    # Решаем
    sol = solve_PDE(sys, g, :LOC)
    @time sol = solve_PDE(sys, g, :LOC)

    # Строим графики
    plot_gif(sol, filename="output/two_dim_example.gif")
    # plot_grid(g, filename="output/two_dim_example_grid.pdf")

end