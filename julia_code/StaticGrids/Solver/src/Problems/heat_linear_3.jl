export solve_PDE

# НЕ ПРОВЕРЕНО!
function solve_PDE(sys::Heat_linear_3, g::UniformGrid)
    # Матрица решение
    u = zeros(g.Nx, g.Nt)

    # Начальное условие
    initial_cond!(u, sys, g)

    # Считаем коэффициенты ур-ния на сетке
    c = zeros(g.Nx, g.Nt - 1)
    k = zeros(g.Nx, g.Nt - 1)
    r = zeros(g.Nx, g.Nt - 1)
    q = zeros(g.Nx, g.Nt - 1)
    f = zeros(g.Nx, g.Nt - 1)
    calculate_koeff!(c, k, r, q, f, sys, g)


    # Простая реализация
    A = zeros(g.Nx - 1)
    B = zeros(g.Nx)
    C = zeros(g.Nx - 1)

    a_1::Float64 = 0
    a_2::Float64 = 0
    b_1::Float64 = 0
    b_2::Float64 = 0
    τ_c::Float64 = 0

    # Цикл по временным слоям
    for j = 1:Nt-1
        τ_c = t_UniformGrid[j] + τ/2


        # Строим тридиагональную матрицу
        B[1] = g.h/(2*g.τ) + sys.β1(τ_c) + sys.α1(τ_c) / h
        C[1] = - sys.α1(τ_c) / g.h
        D[1] = u[1, 1] * (g.h / (2*g.τ)) + sys.μ1(τ_c) + (1/(2*h)) * sys.f(g.x[1], τ_c)
        for i = 2:Nx-1
            a = 0.5 * (k[i-1, j] + k[i, j])
            a_plus = 0.5*(k[i, j] + k[i+1, j])
            R = h * abs(r[i, j]) / (2 * k[i, j])
            ϰ = 1 / (1 + R)
            r_plus = 0.5 * (r[i, j] + abs(r[i, j]))
            r_minus = 0.5 * (r[i, j] - abs(r[i, j]))
            b_plus = r_plus / k[i, j]
            b_minus = r_minus / k[i, j]

            A[i-1] = -a / (g.h^2) * (ϰ - b_minus * g.h)
            C[i-1] = -a_plus / (g.h^2) * (ϰ + b_plus * h)
            B[i] = c[i, j] / g.τ - A[i-1] - C[i-1] - q[i, j]
        end
        B[end] = g.h / (2 * g.τ) + sys.β2(τ_c) + sys.α2(τ_c) / g.h
        A[end] = - sys.α2(τ_c) / h
        D[end] = g.h / (2 * g.τ) * u[end, j] + sys.μ2(τ_c) + (g.h / 2) * f[end, j]

        # Решаем
        M = Tridiagonal(A, B, C)
        u[:, j+1] .= M\D
    end

    return Solution(u, g)
end


function calculate_koeff!(
        c::Vector{Float64},
        k::Vector{Float64},
        r::Vector{Float64},
        q::Vector{Float64},
        f::Vector{Float64},
        sys::Heat_linear_3,
        g::UniformGrid
    )
    τ_c::Float64 = 0
    for j = 1:g.Nt - 1
        τ_c = g.t[j] + g.τ / 2
        for i = 1:g.Nx
            c[i, j] = sys.c(g.x[i], τ_c)
            k[i, j] = sys.k(g.x[i], τ_c)
            r[i, j] = sys.r(g.x[i], τ_c)
            q[i, j] = sys.q(g.x[i], τ_c)
            f[i, j] = sys.f(g.x[i], τ_c)
        end
    end;
end