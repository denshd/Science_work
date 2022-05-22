export solve_PDE

function solve_PDE(sys::Heat2_linear_const_1, g::Grid2, method::Symbol; σ::Tuple{Real, Real}=(1.0, 1.0))
    # Решение
    u = zeros(g.N[1], g.N[2], g.Nt)

    # Учёт начальных/граничных условий

    initial_cond!(u, sys, g)
    boundary_cond_1!(u, sys, g)

    if (method == :LOC)
        loc_scheme!(u, sys, g, σ)
    elseif (method == :sparse)
        sparse_scheme!(u, sys, g)
    end

    return Solution2(u, g)
end


function loc_scheme!(u::Array{Float64, 3}, sys::Heat2_linear_const_1, g::Grid2, σ::Tuple{Real, Real})
    # Временное решение (см. алгоритм)
    w = zeros(g.N[1], g.N[2])

    # Диагонали в одномерных задачах
    A::Tuple{Vector{Float64}, Vector{Float64}} = (
        -(g.τ * σ[1]) / (g.h[1]^2) * ones(g.N[1] - 1),
        -g.τ * σ[2] / g.h[2]^2 * ones(g.N[2] - 1)
    )
    B::Tuple{Vector{Float64}, Vector{Float64}} = (
        (1 + 2 * g.τ * σ[1] / g.h[1]^2) * ones(g.N[1]),
        (1 + 2 * g.τ * σ[2] / g.h[2]^2) * ones(g.N[2])
    )
    C::Tuple{Vector{Float64}, Vector{Float64}} = (
        -(g.τ * σ[1]) / (g.h[1]^2) * ones(g.N[1] - 1),
        -g.τ * σ[2] / g.h[2]^2 * ones(g.N[2] - 1)
    )
    D::Tuple{Vector{Float64}, Vector{Float64}} = (
        Vector{Float64}(undef, g.N[1]),
        Vector{Float64}(undef, g.N[2])
    )
    B[1][1] = B[1][end] = B[2][1] = B[2][end] = 1
    C[1][1] = C[2][1] = A[1][end] = A[2][end] = 0

    # Временные векторы для прогонки
    α::Tuple{Vector{Float64}, Vector{Float64}} = (
        Vector{Float64}(undef, g.N[1]),
        Vector{Float64}(undef, g.N[2])
    )
    β::Tuple{Vector{Float64}, Vector{Float64}} = (
        Vector{Float64}(undef, g.N[1]),
        Vector{Float64}(undef, g.N[2])
    )
    temp_x::Tuple{Vector{Float64}, Vector{Float64}} = (
        Vector{Float64}(undef, g.N[1]),
        Vector{Float64}(undef, g.N[2])
    )

    # Цикл по временным слоям
    for it = 1 : g.Nt - 1
        t_half = g.t[it] + g.τ / 2

        # Сначала решаем задачу вдоль оси x1
        for i2 = 1 : g.N[2]

            # Создаём трёхдиагональную матрицу
            D[1][1] = sys.mu_1(g.x[2][i2], t_half)
            for i1 = 2 : g.N[1] - 1
                D[1][i1] = u[i1, i2, it] + g.τ * (1 - σ[1]) / g.h[2]^2 * (
                    u[i1 + 1, i2, it] - 2 * u[i1, i2, it] + u[i1 - 1, i2, it]
                ) + (g.τ / 2) * sys.f(g.x[1][i1], g.x[2][i2], t_half)
            end
            D[1][end] = sys.mu1(g.x[2][i2], t_half)

            # Получаем временный вектор v_(1)
            tridiagonal_matrix_algorithm!(α[1], β[1], temp_x[1], A[1], B[1], C[1], D[1])
            w[:, i2] .= temp_x[1]
        end

        # Теперь решаем локально-одномерную задачу вдоль оси x2 для каждого x[1][i]:
        for i1 = 1 : g.N[1]

            # Создаём трёхдиагональную матрицу
            D[2][1] = u[i1, 1, it + 1]
            for i2 = 2 : g.N[2] - 1
                D[2][i2] = w[i1, i2] + g.τ * (1 - σ[2]) / g.h[2]^2 * (
                    w[i1, i2 + 1] - 2 * w[i1, i2] + w[i1, i2 - 1]
                ) + g.τ / 2 * sys.f(g.x[1][i1], g.x[2][i2], g.t[it + 1])
            end
            D[2][end] = u[i1, end, it + 1]

            # Получаем решение на временном срезе it + 1 для x[1] = x[1][i]
            tridiagonal_matrix_algorithm!(α[2], β[2], temp_x[2], A[2], B[2], C[2], D[2])
            u[i1, :, it + 1] = temp_x[2]
        end
    end
end


function sparse_scheme!(u::Array{Float64, 3}, sys::Heat2_linear_const_1, g::Grid2)
    # Считаем неоднородность f(x, t) на сетке
    f = Array{Float64, 3}(undef, g.N[1], g.N[2], g.Nt)
    for i = 1:g.N[1]
        for j = 1:g.N[2]
            for k = 1:g.Nt
                f[i, j, k] = sys.f(g.x[1][i], g.x[2][j], g.t[k])
            end
        end
    end

    throw(ErrorException("Ещё не реализовано, соре"))
end
