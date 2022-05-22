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


    function Memory_for_time_step(N::Tuple{Integer, Integer})
        new(
            Tuple(Vector{Cdouble}(undef, N[i] - 1) for i = 1:2),
            Tuple(Vector{Cdouble}(undef, N[i]) for i = 1:2),
            Tuple(Vector{Cdouble}(undef, N[i] - 1) for i = 1:2),
            Tuple(Vector{Cdouble}(undef, N[i]) for i = 1:2),
            Tuple(Vector{Cdouble}(undef, N[i]) for i = 1:2),
            Tuple(Vector{Cdouble}(undef, N[i]) for i = 1:2),
            Tuple(Vector{Cdouble}(undef, N[i]) for i = 1:2)
        )
    end
end


"""
"Блок". Хранит сетку, текущее время, временной шаг блока
"""
struct Block{T <: AbstractRange{Cdouble}}
    spacial_grid::UniformGrid2{T} # Сетка

    block_index::Int
    supblock_index::Int
    supblock_position::Tuple{Tuple{Int, Int}, Tuple{Int, Int}} # Координаты `((x1, x2), (y1, y2))` положения блока в надблоке (в индексных координатах надблока)
    
    u_new::Matrix{Cdouble} # Новое значение
    u_old::Matrix{Cdouble} # Старое значение (просто нельзя вот так взять и обновить новое без старого. А каждый раз выделять память -- затратно)
    w::Matrix{Cdouble} # Промежуточное значение в Локально-одномерной схеме
    w_boundary::Matrix{Cdouble} # Граничные значения промежуточного значения `w` в Локально-одномерной схеме.
    
    memory_for_time_step::Memory_for_time_step # А тут всё то, что нужно для одного time_step!
end



function create_block(
    spacial_grid::UniformGrid2{T},
    block_index::Integer,
    supblock_index::Integer,
    supblock_position::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}}
    ) where T

    Nx = spacial_grid.g[1].N
    Ny = spacial_grid.g[2].N

    return Block{T}(
        spacial_grid,
        block_index,
        supblock_index,
        supblock_position,
        Matrix{Cdouble}(undef, (Nx, Ny)),
        Matrix{Cdouble}(undef, (Nx, Ny)),
        Matrix{Cdouble}(undef, (Nx, Ny)),
        Matrix{Cdouble}(undef, (2, Ny)),
        Memory_for_time_step((Nx, Ny))
    )
end


function create_subblock(
    block::Block{T},
    subblock_index::Integer,
    subblock_position::Tuple{Tuple{Integer, Integer}, Tuple{Integer, Integer}}
    ) where T

    i1 = subblock_position[1][1]
    i2 = subblock_position[1][2]
    j1 = subblock_position[2][1]
    j2 = subblock_position[2][2]

    x = block.spacial_grid.g[1].x
    hx = block.spacial_grid.g[1].h
    y = block.spacial_grid.g[2].x
    hy = block.spacial_grid.g[2].h

    return create_block(
        UniformGrid2(
            ((x[i1], x[i2]), (y[j1], y[j2])), 
            (0.5 * hx, 0.5 * hy)
        ),
        subblock_index,
        block.block_index,
        subblock_position
    )
end


"""
Структура, представляющая собой "уровень"
"""
mutable struct Level{T <: AbstractRange{Cdouble}}
    level_number::Int # номер уровня в иерархии уровней
    sublevel::Union{Nothing, Level{T}} # подуровень
    suplevel::Union{Nothing, Level{T}} # надуровень

    blocks::Vector{Block{T}} # массив блоков сеток UniformGrid2
    M::Int # число блоков сетки UniformGrid2

    t_curr::Cdouble # Текущее время на уровне
    Δt::Cdouble # Временной шаг уровня
end


function initialize_level(blocks::Vector{Block{T}}, t_curr::Cdouble, Δt::Cdouble) where T
    return Level{T}(
        1,
        nothing,
        nothing,
        blocks,
        length(blocks),
        t_curr,
        Δt
    )
end


function add_sublevel!(suplevel::Level{T}, blocks::Vector{Block{T}}) where T
    if has_sublevel(suplevel)
        throw(Exception("level already have sublevel"))
    else
        suplevel.sublevel = Level{T}(
            suplevel.level_number + 1,
            nothing,
            suplevel,
            blocks,
            length(blocks),
            suplevel.t_curr,
            suplevel.Δt * 0.5
        )
    end
    return suplevel.sublevel
end


"""
Функция проверяет, ∃? подуровень уровня `L::Level`
"""
function has_sublevel(L::Level{T}) where T
    if (typeof(L.sublevel) == Nothing)
        return false
    else
        return true
    end
end


"""
Функция проверяет, ∃? надуровень уровня `L::Level`
"""
function has_suplevel(L::Level{T}) where T
    if (typeof(L.suplevel) == Nothing)
        return false
    else
        return true
    end
end


"""
Возвращает `n`-ый подуровень (`0` -- сам уровень) уровня `L`
"""
function get_level(L::Level{T}, n::Int) where T
    if (n == 0)
        return L
    elseif has_sublevel(L)
        return get_level(L.sublevel, n - 1)
    else
        throw(DomainError("$n sublevel doesn't exist."))
    end
end


function get_suplevel(L::Level{T}) where T
    if has_suplevel(L)
        return L.suplevel
    else
        throw(DomainError("suplevel doesn't exist"))
    end
end


function get_sublevel(L::Level{T}) where T
    if has_sublevel(L)
        return L.sublevel
    else
        throw(DomainError("sublevel doesn't exist"))
    end
end


function get_level_refinement_ratio(L::Level{T}) where T
    if has_suplevel(L)
        return 2
    else
        return 1
    end
end


function get_type(L::Level{T}) where T
    return T
end


function get_deepest_level(L::Level{T}) where T
    l = L
    while has_sublevel(L)
        l = L.sublevel
    end
    return l
end


# function create_solution_function(level::Leve{T}) where T
#     curr_level = get_deepest_level(level)
#     function result(x, y, t)
#         while true
#             for block in curr_level.blocks
#                 if x in 
#         end
#     end
# end


"""
Вот здесь хранится ВСЁ.
"""
mutable struct Solution{T <: AbstractRange{Cdouble}}
    levels::Vector{Level{T}} # Вектор нулевых уровней (для каждой временной точки). Нулевой уровень -- он же содержит и подуровни
    time_grid::UniformGrid{T}
end
