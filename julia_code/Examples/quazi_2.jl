using Solver

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

    N = (101, 101)
    Nt = 501

    g = UniformGrid_2((L1, L2), T, N, Nt)

    sys = Heat2_quazilinear_1(k1, k2, f, mu_1, mu1, mu_2, mu2, u0)

    @time sol = solve_PDE(sys, g, :LOC)

    plot_heatmap_gif(sol, filename="Samarski_2.gif")

    plot_gif(sol)
end