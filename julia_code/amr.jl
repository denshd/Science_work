# Попытка самой простой реализации, которая приходит в голову
# using Plots, BenchmarkTools, LinearAlgebra, LaTeXStrings
double = Float64


"""
Двумерная пространственная сетка
"""
struct UniformGrid
    N::Tuple{Int, Int} # Число точек по осям 1 и 2 соотв.
    h::Tuple{double, double} # Шаги по осям 1 и 2 соотв.
    x::Tuple{Vector{double}, Vector{double}} # Точки сеток по осям 1 и 2 соотв.


    """
    Создаёт равномерную сетку на [0, 1]×[0, 1] с числом точек `N`
    """
    function UniformGrid(N::Tuple{Integer, Integer})
        x1 = range(start=0, stop=1, length=N[1])
        x2 = range(start=0, stop=1, length=N[2])
        h = (x1.step, x2.step)

        new(N, h, (x1, x2))
    end


    """
    Создаёт равномерную сетку на [0, 1] × [0, 1] с шагами `h`
    """
    function UniformGrid(h::Tuple{Real, Real})
        x1 = range(start=0, stop=1, step=h[1])
        x2 = range(start=0, stop=1, step=h[2])

        new((x1.len, x2.len), h, (x1, x2))
    end


    """
    Создаёт равномерную сетку на `x_b[1] × x_b[2]` с числом точек `N`
    """
    function UniformGrid(x_b::Tuple{Tuple{Real, Real}, Tuple{Real, Real}}, N::Tuple{Int, Int})
        x1 = range(start=x_b[1][1], stop=x_b[1][2], length=N[1])
        x2 = range(start=x_b[2][1], stop=x_b[2][2], length=N[2])

        new(N, (x1.step, x2.step), (x1, x2))
    end


    """
    Создаёт равномерную сетку на `x_b[1] × x_b[2]` с шагами `h`
    """
    function UniformGrid(x_b::Tuple{Tuple{Real, Real}, Tuple{Real, Real}}, h::Tuple{Real, Real})
        x1 = range(start=x_b[1][1], stop=x_b[1][2], step=h[1])
        x2 = range(start=x_b[2][1], stop=x_b[2][2], step=h[2])

        new((x1.len, x2.len), h, (x1, x2))
    end


    """
    Создаёт равномерную сетку на `[0, L[1]] × [0, L[2]]` с числом точек `N` 
    """
    function UniformGrid(L::Tuple{Real, Real}, N::Tuple{Integer, Integer})
        x1 = range(start=0, stop=L[1], length=N[1])
        x2 = range(start=0, stop=L[2], length=N[2])
        
        new(N, (x2.step, x2.step), (x1, x2))
    end


    """
    Создаёт равномерную сетку на `[0, L[1]] × [0, L[2]]` с шагами `h`
    """
    function UniformGrid(L::Tuple{Real, Real}, h::Tuple{Real, Real})
        x1 = range(start=0, stop=L[1], step=h[1])
        x2 = range(start=0, stop=L[2], step=h[2])

        new((x1.len, x2.len), h, (x1, x2))
    end
end


"""
Двумерная квазилинейная задача теплопроводности
"""
struct HeatProblem{K1<:Function, K2<:Function, Mu_1<:Function, Mu1<:Function, Mu_2<:Function, Mu2<:Function, U<:Function}
    k1::K1
    k2::K2
    mu_1::Mu_1
    mu1::Mu1
    mu_2::Mu_2
    mu2::Mu2
    u0::U
end


"""Простая статическая блочно-структурированная сетка, 2 уровня"""
struct SimpleGrid
    coarse_grid::UniformGrid # Главная, грубая сетка
    fine_grids::Vector{UniformGrid} # Массив мелких сеток

    """
    Простейший конструктор, вручную задающий положение всех `fine_grids`
    
    В алгоритме предполагается, что fine-уровни "корректно" вложены!
    Всякие проверки не реализованы (не до этого ещё)
    """
    function SimpleGrid(coarse::UniformGrid, positions::Vector{<:Tuple{Tuple{Real, Real}, Tuple{Real, Real}}})
        fines = [
            UniformGrid(pos, (coarse.h[1] * 0.5, coarse.h[2] * 0.5)) for pos in positions
        ]
        new(coarse, fines)
    end

end



function plot_grid(g::SimpleGrid, filename="grid.pdf")
    # Построение точечек
    p = scatter(
        [(g.coarse_grid.x[1][i], g.coarse_grid.x[2][j]) for i in 1:g.coarse_grid.N[1] for j in 1:g.coarse_grid.N[2]],
        label="",
        color=:black,
        markersize=3,
        markeralpha=0.5,
        markerstrokecolor=:black
    )
    for fine_grid in g.fine_grids
        scatter!(
            [(fine_grid.x[1][i], fine_grid.x[2][j]) for i in 1:fine_grid.N[1] for j in 1:fine_grid.N[2]],
            label="",
            color=:blue,
            markersize=3,
            markeralpha=0.5,
            markerstrokecolor=:blue
        )
    end

    # Построение линий
    for i = 1:g.coarse_grid.N[1]
        plot!(
            [g.coarse_grid.x[1][i], g.coarse_grid.x[1][i]], [g.coarse_grid.x[2][1], g.coarse_grid.x[2][end]],
            label="",
            color=:black
        )
    end
    for i = 1:g.coarse_grid.N[2]
        plot!(
            [g.coarse_grid.x[1][1], g.coarse_grid.x[1][end]], [g.coarse_grid.x[2][i], g.coarse_grid.x[2][i]],
            label="",
            color=:black
        )
    end

    for fine_grid in g.fine_grids
        for i = 1:fine_grid.N[1]
            plot!(
                [fine_grid.x[1][i], fine_grid.x[1][i]], [fine_grid.x[2][1], fine_grid.x[2][end]],
                label="",
                color=:blue
            )
        end
        for i = 1:fine_grid.N[2]
            plot!(
                [fine_grid.x[1][1], fine_grid.x[1][end]], [fine_grid.x[2][i], fine_grid.x[2][i]],
                label="",
                color=:blue
            )
        end
    end
    savefig(p, filename);
end


"""
Простейший алгоритм решения ур-ния теплопроводности с использованием `SimpleGrid`
"""
# function solve_PDE(sys::Heat_problem, 