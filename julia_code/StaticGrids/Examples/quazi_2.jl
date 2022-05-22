using Solver, Plots

function quazi_2()
    function u(x1, x2, t)
        if (t > (x1 + 2*x2))
            return 0.5 * √(-1 + √(1 + 16 * (t - x1 - 2*x2)))
        else
            return 0
        end
    end

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

    u0(x1, x2) = 0.0
    k1(u) = 4 * u^4
    k2(u) = 0.25 * u^2
    f(u) = 0.0

    L1 = 30
    L2 = 20
    T = 50

    N = (301, 201)
    Nt = 501

    g = UniformGrid_2((L1, L2), T, N, Nt)

    sys = Heat2_quazilinear_1(k1, k2, f, mu_1, mu1, mu_2, mu2, u0)

    @time sol = solve_PDE(sys, g, :LOC)

    # plot_heatmap_gif(sol, filename="output/Samarski_2.gif")


    # plot(sol.g.x[1], sol.u[:, 110, 301])
    # plot!(sol.g.x[1], sol.u[:, 70, 301])
    # p1 = plot!(sol.g.x[1], sol.u[:, 30, 301])
    # savefig(p1, "output/Samarski_2_x.pdf")


    # plot(sol.g.x[2], sol.u[10, :, 301])
    # plot!(sol.g.x[2], sol.u[110, :, 301])
    # p2 = plot!(sol.g.x[2], sol.u[210, :, 301])
    # savefig(p1, "output/Samarski_2_y.pdf")

    # display(plot(p1, p2))

    # plot_gif(sol)


    # display(plot(sol.g.x[1], sol.u[:, 30, 301]))
    return nothing
end