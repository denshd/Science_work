export Problem, Parabolic_problem, Hyperbolic_problem, Elliptic_problem, Heat_linear_const_1, Heat_linear_3, Heat_quazilinear_1

abstract type Problem end
    abstract type Parabolic_problem <: Problem end
    abstract type Hyperbolic_problem <: Problem end
    abstract type Elliptic_problem <: Problem end


    """
    Первая краевая задача для линейного уравнения теплопроводности с постоянными коэффициентам
    """
    struct Heat_linear_const_1{F<:Function, Mu1<:Function, Mu2<:Function, U<:Function} <: Parabolic_problem
        f::F
        mu1::Mu1
        mu2::Mu2 
        u0::U
    end


    """
    Третья краевая задача для линейного уравнения теплопроводности общего вида с переменными коэффициентами
    """
    struct Heat_linear_3{C<:Function,
            K<:Function,
            R<:Function,
            Q<:Function,
            F<:Function,
            U0<:Function,
            A1<:Function,
            A2<:Function,
            B1<:Function,
            B2<:Function,
            Mu1<:Function,
            Mu2<:Function} <: Parabolic_problem
        c::C
        k::K
        r::R
        q::Q
        f::F
        u0::U0
        α1::A1
        α2::A2
        β1::B1
        β2::B2
        μ1::Mu1
        μ2::Mu2
    end

    
    """
    Первая краевая задача для квазилинейного уравнения теплопроводности
    """
    struct Heat_quazilinear_1{K<:Function, F<:Function, Mu1<:Function, Mu2<:Function, U<:Function} <: Parabolic_problem
        k::K
        f::F
        mu1::Mu1
        mu2::Mu2 
        u0::U
    end
#-


"""
Учёт начальных условий
"""
function initial_cond!(u::Matrix{Float64}, sys::Parabolic_problem, g::Grid)
    @. u[:, 1] = sys.u0(g.x)
end


"""
Учёт граничных условий краевой задачи 1-го рода
"""
function boundary_cond_1!(u::Matrix{Float64}, sys::Parabolic_problem, g::Grid)
    @. u[1, :] = sys.mu1(g.t)
    @. u[end, :] = sys.mu2(g.t)
end


include("heat_linear_3.jl")
include("heat_linear_const_1.jl")
include("heat_quazilinear_1.jl")