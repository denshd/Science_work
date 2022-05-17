export Grid, Grid1, Grid2, UniformGrid, GeneralGrid, UniformGrid_2, NGrid, TimeNGrid, UniformNGrid, UniformTime, UniformTimeNGrid


"""
Тип, описывающий "сетки"  
"""
abstract type Grid end

    abstract type Grid1 <: Grid end
    abstract type Grid2 <: Grid end

    """
    `UniformGrid` --- равномерная сетка (x, t)  
    """
    struct UniformGrid <: Grid1
        Nx::Int64
        Nt::Int64
        h::Float64
        τ::Float64
        x::Vector{Float64}
        t::Vector{Float64}

        UniformGrid(Nx::Integer, Nt::Integer) = begin
            x = range(start=0, stop=1, length=Nx)
            t = range(start=0, stop=1, length=Nt)
            new(Nx, Nt, x.step, t.step, x, t)
        end

        UniformGrid(h::Real, τ::Real) = begin
            x = range(start=0, stop=1, step=h)
            t = range(start=0, stop=1, step=τ)
            new(x.len, t.len, h, τ, x, t)
        end

        UniformGrid((x1, x2)::Tuple{Real, Real}, (t1, t2)::Tuple{Real, Real}, Nx::Integer, Nt::Integer) = begin
            x = range(start=x1, stop=x2, length=Nx)
            t = range(start=t1, stop=t2, length=Nt)
            new(Nx, Nt, x.step, t.step, x, t)
        end

        UniformGrid((x1, x2)::Tuple{Real, Real}, (t1, t2)::Tuple{Real, Real}, h::Real, τ::Real) = begin
            x = range(start=x1, stop=x2, step=h)
            t = range(start=t1, stop=t2, step=τ)
            new(x.len, t.len, h, τ, x, t)
        end

        UniformGrid(L::Real, T::Real, Nx::Integer, Nt::Integer) = begin
            x = range(start=0, stop=L, length=Nx)
            t = range(start=0, stop=T, length=Nt)
            new(Nx, Nt, x.step, t.step, x, t)
        end

        UniformGrid(L::Real, T::Real, h::Real, τ::Real) = begin
            x = range(start=0, stop=L, step=h)
            t = range(start=0, stop=T, step=τ)
            new(x.len, t.len, h, τ, x, t)
        end
    end


    """
    `GeneralGrid` --- произвольная неравномерная сетка (x, t)
    """
    struct GeneralGrid <: Grid1
        Nx::Int64
        Nt::Int64
        x::Vector{Float64}
        t::Vector{Float64}
        h::Vector{Float64}
        τ::Vector{Float64}

        function GeneralGrid(x::Vector{<:Real}, t::Vector{<:Real})
            Nx = length(x)
            Nt = length(t)
            h = [x[i+1] - x[i] for i = 1:Nx-1]
            τ = [t[i+1] - t[i] for i = 1:Nt-1]

            new(Nx, Nt, x, t, h, τ)
        end
    end


    struct UniformGrid_2{T <: AbstractRange{Float64}} <: Grid2
        N::Tuple{Int, Int}
        Nt::Int
        h::Tuple{Float64, Float64}
        τ::Float64
        x::Tuple{T, T}
        t::T

        function UniformGrid_2(N::Tuple{Integer, Integer}, Nt::Integer)
            x1 = range(start=0, stop=1, length=N[1])
            x2 = range(start=0, stop=1, length=N[2])
            t = range(start=0, stop=1, length=Nt)

            h = (x1.step, x2.step)
            τ = t.step

            new{typeof(x1)}(N, Nt, h, τ, (x1, x2), t)
        end

        function UniformGrid_2(h::Tuple{Real, Real}, τ::Real)
            x1 = range(start=0, stop=1, step=h[1])
            x2 = range(start=0, stop=1, step=h[2])
            t = range(start=0, stop=1, step=τ)

            new{typeof(x1)}((x1.len, x2.len), t.len, h, τ, (x1, x2), t)
        end

        function UniformGrid_2(x_b::Tuple{Tuple{Real, Real}, Tuple{Real, Real}}, t_b::Tuple{Real, Real}, N::Tuple{Int64, Int64}, Nt::Integer)
            x1 = range(start=x_b[1][1], stop=x_b[1][2], length=N[1])
            x2 = range(start=x_b[2][1], stop=x_b[2][2], length=N[2])
            t = range(start=t_b[1], stop=t_b[2], length=Nt)

            new{typeof(x1)}(N, Nt, (x1.step, x2.step), t.step, (x1, x2), t)
        end

        # TODO: fix the problem (N isn't define)
        function UniformGrid_2(x_b::Tuple{Tuple{Real, Real}, Tuple{Real, Real}}, t_b::Tuple{Real, Real}, h::Tuple{Real, Real}, τ::Real)
            x1 = range(start=x_b[1][1], stop=x_b[1][2], step=h[1])
            x2 = range(start=x_b[2][1], stop=x_b[2][2], step=h[2])
            t = range(start=t_b[1], stop=t_b[2], step=τ)

            new{typeof(x1)}(N, Nt, h, τ, (x1, x2), t)
        end

        function UniformGrid_2(L::Tuple{Real, Real}, T::Real, N::Tuple{Integer, Integer}, Nt::Integer)
            x1 = range(start=0, stop=L[1], length=N[1])
            x2 = range(start=0, stop=L[2], length=N[2])
            t = range(start=0, stop=T, length=Nt)
            
            new{typeof(x1)}(N, Nt, (x2.step, x2.step), t.step, (x1, x2), t)
        end

        function UniformGrid_2(L::Tuple{Real, Real}, T::Real, h::Tuple{Real, Real}, τ::Real)
            x1 = range(start=0, stop=L[1], step=h[1])
            x2 = range(start=0, stop=L[2], step=h[2])
            t = range(start=0, stop=T, step=τ)

            new{typeof(x1)}((x1.len, x2.len), t.len, h, τ, (x1, x2), t)
        end
    end
#-


"""
Абстрактный тип сетка (имеется ввиду пространственная сетка)
"""
abstract type NGrid end


"""
Абстрактный тип пространственно-временная сетка
"""
abstract type TimeNGrid end


"""
Равномерная пространственная сетка
"""
struct UniformNGrid <: NGrid

    dim::Int
    N::Vector{Int}
    h::Vector{Float64}
    x::Vector{Vector{Float64}}

    function UniformNGrid(dim::Integer, N::Vector{<:Integer})
        if (length(N) != dim)
            throw(DimensionMismatch("dim != length(N)"))
        end

        x::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, dim)
        h::Vector{Float64} = Vector{Float64}(undef, dim)
        for i = 1:dim
            x[i] = collect(range(start=0, stop=1, length=N[i]))
            h[i] = x[i][2] - x[i][1]
        end

        new(dim, N, h, x)
    end


    function UniformNGrid(dim::Integer, h::Vector{<:Real})
        if (length(h) != dim)
            throw(DimensionMismatch("dim != length(h)"))
        end

        x::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, dim)
        N::Vector{Int} = Vector{Int}(undef, dim)

        for i = 1:dim
            x[i] = collect(range(start=0, stop=1, step=h[i]))
            N[i] = length(x[i])
        end

        new(dim, N, h, x)
    end


    function UniformNGrid(x_b::Vector{Tuple{Real, Real}}, N::Vector{<:Integer})
        x::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, dim)
        h::Vector{Float64} = Vector{Float64}(undef, dim)

        for i = 1:dim
            x[i] = collect(range(start=x_b[i][1], stop=x_b[i][2], length=N[i]))
            h[i] = x[i][2] - x[i][1]
        end

        new(dim, N, h, x)
    end


    function UniformNGrid(x_b::Vector{Tuple{Real, Real}}, h::Vector{Float64})
        x::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, dim)
        N::Vector{Int} = Vector{Int}(undef, dim)

        for i = 1:dim
            x[i] = collect(range(start=x_b[i][1], stop=x_b[i][2], step=h[i]))
            N[i] = length(x[i])
        end

        new(dim, N, h, x)
    end


    function UniformNGrid(L::Vector{<:Real}, N::Vector{<:Integer})
        x::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, dim)
        h::Vector{Float64} = Vector{Float64}(undef, dim)

        for i = 1:dim
            x[i] = collect(range(start=0, stop=L[i], length=N[i]))
            h[i] = x[i][2] - x[i][1]
        end

        new(dim, N, h, x)
    end


    function UniformNGrid(L::Vector{<:Real}, h::Vector{<:Real})
        x::Vector{Vector{Float64}} = Vector{Vector{Float64}}(undef, dim)
        N::Vector{Int} = Vector{Int}(undef, dim)

        for i = 1:dim
            x[i] = collect(range(start=0, stop=L[i], step=h[i]))
            N[i] = length(x[i])
        end

        new(dim, N, h, x)
    end
end


"""
Равномерная сетка по времени
"""
struct UniformTime
    N::Int
    τ::Float64
    t::Vector{Float64}


    function UniformTime(N::Integer)
        t = range(start=0, stop=1, length=N)
        τ = t.step

        new(N, τ, t)
    end


    function UniformTime(τ::Float64)
        t = range(start=0, stop=1, step=τ)
        N = t.len

        new(N, τ, t)
    end


    function UniformTime(t_b::Tuple{Real, Real}, N::Integer)
        t = range(start=t_b[1], stop=t_b[2], length=N)
        τ = t.step

        new(N, τ, t)
    end


    function UniformTime(t_b::Tuple{Real, Real}, τ::Float64)
        t = range(start=_b[1], stop=t_b[2], step=τ)
        N = t.len

        new(N, τ, t)
    end


    function UniformTime(T::Real, N::Integer)
        t = range(start=0, stop=T, length=N)
        τ = t.step

        new(N, τ, t)
    end


    function UniformTime(T::Real, τ::Real)
        t = range(start=0, stop=T, step=τ)
        N = t.len

        new(N, τ, t)
    end
end


"""
Равномерная пространственно-временная сетка
"""
struct UniformTimeNGrid <: TimeNGrid
    g::UniformNGrid
    t::UniformTime
end


    # """
    # `Agrid` --- адаптивная блочная сетка (x, t)
    # """
    # struct AGrid <: Grid
    #     coarse::UniformGrid
    #     subdivision_level::Int64
    #     fine_bounds::Vector{Tuple{Int64, Int64}}
    # end
