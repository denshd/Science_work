include("schemes.jl")

using Plots, GLMakie, LaTeXStrings, LinearAlgebra, BenchmarkTools, Printf
gr()

# Первая функция
function test_1()
    f(x, t) = 0
    mu1(t) = 0
    mu2(t) = 0
    u0(x) = x * sin(3π * x)

    T = 0.05
    Nx = 501
    Nt = 25_001

    @benchmark u = explicit_scheme($f, $u0, $mu1, $mu2, $T, $Nx, $Nt)

    # Plots.plot(range(0, 1, length = Nx), range(0, T, length = Nt), u', st=:wireframe, 
    # xlabel=L"x",
    # ylabel=L"t",
    # zlabel=L"u(x, t)")
end


function test_2()
    k(u) = 1
    f(x, t) = 0
    mu1(t) = 0
    mu2(t) = 0
    u0(x) = x * sin(3π * x)

    T = 0.05
    Nx = 11
    Nt = 501

    u = implicit_scheme(k, f, u0, mu1, mu2, T, Nx, Nt)

    Plots.plot(range(0, 1, length = Nx), range(0, T, length = Nt), u')
end