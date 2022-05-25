include("schemes.jl")

using Plots, GLMakie, LaTeXStrings, LinearAlgebra, BenchmarkTools, Printf
pgfplotsx()


function analytical_solution(x, t)
    res = 0.0
    nfinity::Int = 200
    for m in 1:nfinity
        res += m * ℯ^(-(2π * m)^2 * t) * sin(2π * m * x) / (9 - 4m^2)^2
    end
    res *= -48 / π^2
    res += 0.5 * ℯ^(-(3π)^2 * t) * sin(3π * x)
    return res
end


function numerical_surface()
    f(x, t) = 0
    mu1(t) = 0
    mu2(t) = 0
    u0(x) = x * sin(3π * x)

    T = 0.05
    Nx = 501
    Nt = 25_001
    x = range(0, 1, length = Nx)
    t = range(0, T, length = Nt)

    println("Solving...")
    u = explicit_scheme(f, u0, mu1, mu2, T, Nx, Nt)
    println("done!")

    println("Plotting...")
    number_of_t_dots = 50
    number_of_x_dots = 50
    stepik_t = Nt ÷ number_of_t_dots
    stepik_x = Nx ÷ number_of_x_dots
    p = Plots.plot(
        x[1:stepik_x:end], t[1:stepik_t:end], u[1:stepik_x:end, 1:stepik_t:end]',
        label="",
        titlefont=(14, "Computer Modern"),
        tickfont=(10, "Computer Modern"),
        guidefont=(14, "Computer Modern"),
        xlabel=L"x",
        ylabel=L"t",
        zlabel=L"u_h(x, t)",
        st = :surface,
        camera=(210, 30),
        cbar=false
    )
    println("done!")

    println("Saving...")
    Plots.savefig(p, joinpath("output/", "problem_1_explicit_surface.pdf"))
    println("done!")

    return nothing
end


function analytical_surface()
    f(x, t) = 0
    mu1(t) = 0
    mu2(t) = 0
    u0(x) = x * sin(3π * x)

    T = 0.05
    Nx = 101
    Nt = 101
    x = range(0, 1, length = Nx)
    t = range(0, T, length = Nt)

    println("Solving...")
    u = zeros(Nx, Nt)
    for i in 1:Nx
        for j in 1:Nt
            u[i, j] = analytical_solution(x[i], t[j])
        end
    end
    println("done!")

    println("Plotting...")
    number_of_t_dots = 50
    number_of_x_dots = 50
    stepik_t = Nt ÷ number_of_t_dots
    stepik_x = Nx ÷ number_of_x_dots
    p = Plots.plot(
        x[1:stepik_x:end], t[1:stepik_t:end], u[1:stepik_x:end, 1:stepik_t:end]',
        label="",
        titlefont=(14, "Computer Modern"),
        tickfont=(10, "Computer Modern"),
        guidefont=(14, "Computer Modern"),
        xlabel=L"x",
        ylabel=L"t",
        zlabel=L"u(x, t)",
        st = :surface,
        camera=(210, 30),
        cbar=false
    )
    println("done!")

    println("Saving...")
    Plots.savefig(p, joinpath("output/", "problem_1_analytic_surface.pdf"))
    println("done!")

    return nothing
end


function error_surface()
    f(x, t) = 0
    mu1(t) = 0
    mu2(t) = 0
    u0(x) = x * sin(3π * x)

    T = 0.05
    Nx = 501
    Nt = 25_001
    x = range(0, 1, length = Nx)
    t = range(0, T, length = Nt)

    println("Solving numerically...")
    u_h = explicit_scheme(f, u0, mu1, mu2, T, Nx, Nt)
    println("done!")

    println("Solving analytically...")
    u = zeros(Nx, Nt)
    for i in 1:Nx
        for j in 1:Nt
            u[i, j] = analytical_solution(x[i], t[j])
        end
    end
    println("done!")
    
    println("Calculating local error...")
    err = zeros(Nx, Nt)
    for i in 1:Nx
        for j in 1:Nt
            err[i, j] = abs(u[i, j] - u_h[i, j])
        end
    end
    println("done!")

    println("Plotting...")
    number_of_t_dots = 50
    number_of_x_dots = 50
    stepik_t = Nt ÷ number_of_t_dots
    stepik_x = Nx ÷ number_of_x_dots
    p = Plots.plot(
        x[1:stepik_x:end], t[1:stepik_t:end], err[1:stepik_x:end, 1:stepik_t:end]',
        label="",
        titlefont=(14, "Computer Modern"),
        tickfont=(10, "Computer Modern"),
        guidefont=(14, "Computer Modern"),
        xlabel=L"x",
        ylabel=L"t",
        zlabel=L"|u(x, t) - u_h(x, t)|",
        st = :surface,
        camera=(200, 30),
        cbar=false
    )
    println("done!")

    println("Saving...")
    Plots.savefig(p, joinpath("output/", "problem_1_explicit_error_surface.pdf"))
    println("done!")

    return nothing
end


function time_dependence()
    f(x, t) = 0
    mu1(t) = 0
    mu2(t) = 0
    u0(x) = x * sin(3π * x)

    T = 0.05
    Nt = 25_001
    Nx = 11:50:501

    times = Vector{Cdouble}(undef, length(Nx))
    for (i, N) in enumerate(Nx)
        println(@sprintf("%.3f", i / length(Nx)))
        times[i] = mean((@benchmark explicit_scheme($f, $u0, $mu1, $mu2, $T, $N, $Nt)).times) * 10^(-6)
    end


    p = Plots.plot(
        Nx, times,
        label="",
        xlabel=L"N_x",
        ylabel=L"$t(N_x)$ [ms]",
        color=:black
    )

    savefig(p, joinpath("output/", "problem_1_explicit_time.pdf"))
end