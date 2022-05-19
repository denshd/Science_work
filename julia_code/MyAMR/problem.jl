"""
Задача Дирихле для двумерного квазилинейного уравнения теплопроводности
"""
struct HeatProblem2{K1<:Function, K2<:Function, F<:Function, Mu_1<:Function, Mu1<:Function, Mu_2<:Function, Mu2<:Function, U<:Function}
    k1::K1
    k2::K2
    f::F
    mu_1::Mu_1
    mu1::Mu1
    mu_2::Mu_2
    mu2::Mu2
    u0::U
end