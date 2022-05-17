println("Including amr.jl...")
include("amr.jl")
println("done!")

println("Including Plots.jl...")
using Plots
println("done!")


"""Граничное условие при `x1 = 0`"""
function mu_1(x2, t)
    if (t > (2*x2))
        return 0.5 * √(-1 + √(1 + 16 * (t - 2*x2)))
    else
        return 0
    end
end


"""Граничное условие при `x1 = L[1]`"""
function mu1(x2, t)
    if (t > (30 + 2*x2))
        return 0.5 * √(-1 + √(1 + 16 * (t - 30 - 2*x2)))
    else
        return 0
    end
end


"""Граничное условие при `x2 = 0`"""
function mu_2(x1, t)
    if (t > x1)
        return 0.5 * √(-1 + √(1 + 16 * (t - x1)))
    else
        return 0
    end
end


"""Граничное условие при `x2 = L[2]`"""
function mu2(x1, t)
    if (t > (x1 + 40))
        return 0.5 * √(-1 + √(1 + 16 * (t - x1 - 40)))
    else
        return 0
    end
end


"Начальное условие"
u0(x1, x2) = 0.0

"Коэффициент `k1` уравнения"
k1(u) = 4 * u^4

"Коэффициент `k2` уравнения"
k2(u) = 0.25 * u^2


# Задача
prob = HeatProblem(k1, k2, mu_1, mu1, mu_2, mu2, u0)


# Параметры пространственной и временной сеток:

L1 = 30
L2 = 20
T = 50

N = (31, 21)
Nt = 51

# Попробуем задать простейшую блочно-структурированную статическую сетку
g = SimpleGrid(
    UniformGrid((L1, L2), N),
    [
        ((5, 15), (2, 12)),
        ((15, 18), (8, 18)),
        ((18, 29), (10, 12))
    ]
)
plot_grid(g)


