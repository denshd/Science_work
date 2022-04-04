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
struct Solution{G<:Grid}
    u::Matrix{Float64}
    g::G
end

include("Plotter.jl")



end