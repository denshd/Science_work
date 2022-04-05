export Problem, Parabolic_problem, Hyperbolic_problem, Elliptic_problem, Heat_linear_const_1, Heat_linear_3, Heat_quazilinear_1, Heat2_linear_const_1, Heat2_linear_1

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



    struct Heat2_linear_const_1{F<:Function, Mu_1<:Function, Mu1<:Function, Mu_2<:Function, Mu2<:Function, U<:Function} <: Parabolic_problem
        f::F
        mu_1::Mu_1
        mu1::Mu1
        mu_2::Mu_2
        mu2::Mu2
        u0::U
    end

    struct Heat2_linear_1{K1<:Function, K2<:Function, F<:Function,
        Chi_minus_1<:Function, Beta_minus_1<:Function, Mu_minus_1<:Function, Chi_plus_1<:Function, Beta_plus_1<:Function, Mu_plus_1<:Function,
        Chi_minus_2<:Function, Beta_minus_2<:Function, Mu_minus_2<:Function, Chi_plus_2<:Function, Beta_plus_2<:Function, Mu_plus_2<:Function,
        U<:Function} <: Parabolic_problem

        k1::K1
        k2::K2
        f::F

        chi_minus_1::Chi_minus_1
        beta_minus_1::Beta_minus_1
        mu_minus_1::Mu_minus_1
        chi_plus_1::Chi_plus_1
        beta_plus_1::Beta_plus_1
        mu_plus_1::Mu_plus_1
        
        chi_minus_2::Chi_minus_2
        beta_minus_2::Beta_minus_2
        mu_minus_2::Mu_minus_2
        chi_plus_2::Chi_plus_2
        beta_plus_2::Beta_plus_2
        mu_plus_2::Mu_plus_2

        u0::U
    end

#-


"""
Учёт начальных условий
"""
function initial_cond!(u::Matrix{Float64}, sys::Parabolic_problem, g::Grid1)
    @. u[:, 1] = sys.u0(g.x)
end


"""
Учёт граничных условий краевой задачи 1-го рода
"""
function boundary_cond_1!(u::Matrix{Float64}, sys::Parabolic_problem, g::Grid1)
    @. u[1, :] = sys.mu1(g.t)
    @. u[end, :] = sys.mu2(g.t)
end


function initial_cond!(u::Array{Float64, 3}, sys::Parabolic_problem, g::Grid2)
    for i = 1:g.N[1]
        for j = 1:g.N[2]
            u[i, j, 1] = sys.u0(g.x[1][i], g.x[2][j])
        end
    end
end


function boundary_cond_1!(u::Array{Float64, 3}, sys::Parabolic_problem, g::Grid2)
    for j = 1:g.N[2]
        for k = 1:g.Nt
            u[1, j, k] = sys.mu_1(g.x[2][j], g.t[k])
            u[end, j, k] = sys.mu1(g.x[2][j], g.t[k])
        end
    end

    for i = 1:g.N[1]
        for k = 1:g.Nt
            u[i, 1, k] = sys.mu_2(g.x[1][i], g.t[k])
            u[i, end, k] = sys.mu2(g.x[1][i], g.t[k])
        end
    end
end


include("heat_linear_3.jl")
include("heat_linear_const_1.jl")
include("heat_quazilinear_1.jl")
include("heat2_linear_const_1.jl")