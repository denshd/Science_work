export Grid, UniformGrid, GeneralGrid


"""
Тип, описывающий "сетки"  
"""
abstract type Grid end

    """
    `UniformGrid` --- равномерная сетка (x, t)  
    """
    struct UniformGrid <: Grid
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
    struct GeneralGrid <: Grid
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


    # """
    # `Agrid` --- адаптивная блочная сетка (x, t)
    # """
    # struct AGrid <: Grid
    #     coarse::UniformGrid
    #     subdivision_level::Int64
    #     fine_bounds::Vector{Tuple{Int64, Int64}}
    # end
