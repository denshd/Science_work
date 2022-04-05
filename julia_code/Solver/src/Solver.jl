"""
Основной модуль "решатель"
"""
module Solver

export Solution


include("tridiagonal_matrix_algorithm.jl")
include("Grids.jl")
include("Problems/Problems.jl")

"""
Структура для хранения решения и сетки, на котором решение получено
"""
struct Solution{G<:Grid1}
    u::Matrix{Float64}
    g::G
end

struct Solution2{G<:Grid2}
    u::Array{Float64, 3}
    g::G
end

include("Plotter.jl")



end