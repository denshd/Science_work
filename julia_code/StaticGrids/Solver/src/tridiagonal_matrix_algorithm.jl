"""
Метод прогонки для решения трёхдиагональных систем линейных уравнений.
"""
function tridiagonal_matrix_algorithm(a::Vector{<:Number}, B::Vector{<:Number}, c::Vector{<:Number}, D::Vector{<:Number})::Vector{Float64}

    # размер матрицы
    n = length(B)

    # удобная проверка
    if (length(a) != length(B))
        A = vcat([0], copy(a))
    else
        A = copy(a)
    end
    if (length(c) != length(B))
        C = vcat(copy(c), [0])
    else
        C = copy(c)
    end

    # вспомогательный вектор y и коэффициенты алгоритма прогонки α и β, и столбец решение x
    # y = Vector{Float64}(undef, n)
    α = Vector{Float64}(undef, n)
    β = Vector{Float64}(undef, n)
    x = Vector{Float64}(undef, n)

    # Алгоритм прогонки:

    # прямой ход
    α[1] = - C[1] / B[1]
    β[1] = D[1] / B[1]
    for i=2:n-1
        α[i] = - C[i] / (B[i] + A[i] * α[i-1])
        β[i] = (D[i] - A[i] * β[i-1]) / (B[i] + A[i] * α[i-1])
    end
    β[n] = (D[n] - A[n] * β[n-1]) / (B[n] + A[n] * α[n-1])

    # обратный ход
    x[n] = β[n]
    for i = n-1:-1:1
        x[i] = α[i]*x[i+1] + β[i]
    end

    return x
end;


"""
Оптимальный по памяти алгоритм прогонки в случае, если его нужно запускать много раз
"""
function tridiagonal_matrix_algorithm!(
        α::Vector{<:T}, β::Vector{<:T}, x::Vector{<:T},
        A::Vector{<:V}, B::Vector{<:V}, C::Vector{<:V}, D::Vector{<:V}
    ) where {T <: Number, V <: Number}

    α[1] = - C[1] / B[1]
    β[1] = D[1] / B[1]
    @inbounds for i=2:(length(B)-1)
        α[i] = - C[i] / (B[i] + A[i-1] * α[i-1])
        β[i] = (D[i] - A[i-1] * β[i-1]) / (B[i] + A[i-1] * α[i-1])
    end
    β[end] = (D[end] - A[end] * β[end-1]) / (B[end] + A[end] * α[end-1])

    # обратный ход
    x[end] = β[end]
    for i = length(B)-1:-1:1
        x[i] = α[i]*x[i+1] + β[i]
    end;
end