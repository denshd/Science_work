"""
Одномерная равномерная сетка из точечек
"""
struct UniformGrid{T <: AbstractRange{Cdouble}}
    N::Int
    h::Cdouble
    x::T


    function UniformGrid(N::Integer)
        x = range(start=0, stop=1, length=N)

        new{typeof(x)}(N, step(x), x)
    end


    function UniformGrid(h::Real)
        x = range(start=0, stop=1, step=h)

        new{typeof(x)}(length(x), h, x)
    end


    function UniformGrid(L::Real, N::Integer)
        x = range(start=0, stop=L, length=N)

        new{typeof(x)}(N, step(x), x)
    end


    function UniformGrid(L::Real, h::Real)
        x = range(start=0, stop=L, step=h)

        new{typeof(x)}(length(x), h, x)
    end


    function UniformGrid(xb::Tuple{Real, Real}, N::Integer)
        x = range(start=xb[1], stop=xb[2], length=N)

        new{typeof(x)}(N, step(x), x)
    end


    function UniformGrid(xb::Tuple{Real, Real}, h::Real)
        x = range(start=xb[1], stop=xb[2], step=h)

        new{typeof(x)}(length(x), h, x)
    end
end


# В таком формализме TimeGrid~--- просто и есть UniformGrid


"""
Двумерная равномерная сетка из точек
"""
struct UniformGrid2{T <: AbstractRange{Cdouble}}
    g::Tuple{UniformGrid{T}, UniformGrid{T}}


    function UniformGrid2(N::Tuple{Integer, Integer})
        g = Tuple(UniformGrid(N[i]) for i = 1:2)
        new{typeof(g[1].x)}(g)
    end


    function UniformGrid2(h::Tuple{Real, Real})
        g = Tuple(UniformGrid(h[i]) for i = 1:2)
        new{typeof(g[1].x)}(g)
    end


    function UniformGrid2(L::Tuple{Real, Real}, N::Tuple{Integer, Integer})
        g = Tuple(UniformGrid(L[i], N[i]) for i = 1:2)
        new{typeof(g[1].x)}(g)
    end


    function UniformGrid2(L::Tuple{Real, Real}, h::Tuple{Real, Real})
        g = Tuple(UniformGrid(L[i], h[i]) for i = 1:2)
        new{typeof(g[1].x)}(g)
    end


    function UniformGrid2(xb::Tuple{Tuple{Real, Real}, Tuple{Real, Real}}, N::Tuple{Integer, Integer})
        g = Tuple(UniformGrid(xb[i], N[i]) for i = 1:2)
        new{typeof(g[1].x)}(g)
    end


    function UniformGrid2(xb::Tuple{Tuple{Real, Real}, Tuple{Real, Real}}, h::Tuple{Real, Real})
        g = Tuple(UniformGrid(xb[i], h[i]) for i = 1:2)
        new{typeof(g[1].x)}(g)
    end
end


"""
Возвращает этот пресловутый тип T
"""
function get_range_type(::UniformGrid2{T}) where T
    return T
end


struct Memory_for_time_step
    A::Tuple{Vector{Cdouble}, Vector{Cdouble}}
    B::Tuple{Vector{Cdouble}, Vector{Cdouble}}
    C::Tuple{Vector{Cdouble}, Vector{Cdouble}}
    D::Tuple{Vector{Cdouble}, Vector{Cdouble}}

    α::Tuple{Vector{Cdouble}, Vector{Cdouble}}
    β::Tuple{Vector{Cdouble}, Vector{Cdouble}}
    temp_x::Tuple{Vector{Cdouble}, Vector{Cdouble}}
end


"""
"Блок". Хранит сетку, текущее время, временной шаг блока
"""
struct Block{T <: AbstractRange{Cdouble}}
    spacial_grid::UniformGrid2{T} # Сетка
    u_new::Matrix{Cdouble} # Новое значение
    u_old::Matrix{Cdouble} # Старое значение (просто нельзя вот так взять и обновить новое без старого. А каждый раз выделять память -- затратно)
    w::Matrix{Cdouble}
    
    """А тут всё то, что нужно для одного time_step!"""
    memory_for_time_step::Memory_for_time_step
end


"""
Структура, представляющая собой "уровень"
"""
mutable struct Level{T <: AbstractRange{Cdouble}}
    level_number::Int # номер уровня в иерархии уровней
    sublevel::Union{Nothing, Level{T}}

    blocks::Vector{Block{T}} # массив блоков сеток UniformGrid2
    M::Int # число блоков сетки UniformGrid2

    t_curr::Cdouble # Текущее время на уровне
    Δt::Cdouble # Временной шаг уровня


    # function Level(spacial_grids::Vector{UniformGrid2{T}}, time_grid::UniformGrid{T}, level_number::Integer) where T <: AbstractRange{Cdouble}
    #     new{get_range_type(spacial_grid[1])}(
    #         spacial_grids,
    #         length(spacial_grids),
    #         time_grid,
    #         Matrix{Cdouble}(undef, ())
    #     )
    # end
end


"""
Функция проверяет, ∃? подуровень уровня `L::Level`
"""
function has_sublevel(L::Level)
    if (typeof(L.sublevel) == Nothing)
        return true
    else
        return false
    end
end


"""
Возвращает `n`-ый подуровень (`0` -- сам уровень) уровня `L`
"""
function get_level(L::Level, n::Int)
    if (n == 0)
        return L
    elseif has_sublevel(L)
        return get_level(L.sublevel, n - 1)
    else
        throw(DomainError("$n sublevel doesn't exist."))
    end
end


"""
Вот здесь хранится ВСЁ.
"""
mutable struct Solution{T <: AbstractRange{Cdouble}}
    levels::Vector{Level{T}} # Вектор нулевых уровней (для каждой временной точки). Нулевой уровень -- он же содержит и подуровни
    time_grid::UniformGrid{T}
end