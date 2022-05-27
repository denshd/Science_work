include("schemes.jl")

using Plots, GLMakie, LaTeXStrings, LinearAlgebra, BenchmarkTools, Printf
pgfplotsx()


function time_dependence()
    k(u) = 1
    f(x, t) = 0
    mu1(t) = 0
    mu2(t) = 0
    u0(x) = x * sin(3π * x)

    T = 0.05
    Nt = 501
    Nx = 11:50:501

    times = Vector{Cdouble}(undef, length(Nx))
    for (i, N) in enumerate(Nx)
        println(@sprintf("%.3f", i / length(Nx)))
        times[i] = mean((@benchmark implicit_scheme($k, $f, $u0, $mu1, $mu2, $T, $N, $Nt)).times) * 10^(-6)
    end


    p = Plots.plot(
        Nx, times,
        label="",
        xlabel=L"N_x",
        ylabel=L"$t(N_x)$ [ms]",
        color=:black,
        titlefont=(14, "Computer Modern"),
        tickfont=(10, "Computer Modern"),
        guidefont=(14, "Computer Modern"),
    )

    savefig(p, joinpath("output/", "problem_1_implicit_time.pdf"))
end


function numerical_surface()
    k(u) = 0.5 * u^2
    f(x, t) = 0
    mu1(t) = 10 * √(t + 0.1)
    mu2(t) = 0
    u0(x) = (x < 0.5) ? (2 * √(5 * (0.5 - x))) : 0

    T = 0.1
    Nx = 51
    Nt = 501
    x = range(0, 1, length = Nx)
    t = range(0, T, length = Nt)

    println("Solving...")
    u = implicit_scheme(k, f, u0, mu1, mu2, T, Nx, Nt)
    println("done!")

    println("Plotting surface...")
    p = Plots.plot(
        x, t .+ 0.1, u',
        label="",
        titlefont=(14, "Computer Modern"),
        tickfont=(10, "Computer Modern"),
        guidefont=(14, "Computer Modern"),
        xlabel=L"x",
        ylabel=L"t",
        zlabel=L"u_h(x, t)",
        st = :surface,
        camera=(-15, 20),
        colormap_name = "viridis",
        cbar=false
    )
    println("done!")

    println("Saving surface...")
    Plots.savefig(p, joinpath("output/", "problem_2_implicit_surface.pdf"))
    println("done!")


    # println("Plotting like Samarski...")
    # index = [1 101 301 501]
    # xlimits = (0, 1)
    # ylimits = (0, 7)
    # ps = Plots.plot()
    # for i ∈ index
    #     Plots.plot!(ps, x, u[:, i], xlims=xlimits, ylims=ylimits, label=latexstring("t = $(t[i])"))
    #     Plots.scatter!(p, x, u[:, i], label="")
    # end
    # Plots.plot!(
    #     titlefont=(14, "Computer Modern"),
    #     tickfont=(10, "Computer Modern"),
    #     guidefont=(14, "Computer Modern"),
    #     xlabel=L"x",
    #     ylabel=L"u_h(x, t)"
    # )
    # println("done!")

    # println("Saving Samarski...")
    # Plots.savefig(ps, joinpath("output/", "problem_2_implicit_Samarski.pdf"))


    return nothing
end


function numerical_Samarski()
    k(u) = 0.5 * u^2
    f(x, t) = 0
    mu1(t) = 10 * √(t + 0.1)
    mu2(t) = 0
    u0(x) = (x < 0.5) ? (2 * √(5 * (0.5 - x))) : 0

    T = 0.1
    Nx = 51
    Nt = 501
    x = range(0, 1, length = Nx)
    t = range(0, T, length = Nt)

    println("Solving...")
    u = implicit_scheme(k, f, u0, mu1, mu2, T, Nx, Nt)
    println("done!")

    println("Plotting like Samarski...")
    index = [1 101 301 501]
    xlimits = (0, 1)
    ylimits = (0, 5)
    ps = Plots.plot()
    for i ∈ index
        Plots.plot!(
            ps, x, u[:, i], xlims=xlimits, ylims=ylimits, label=latexstring("t = $(t[i])"),
            titlefont=(14, "Computer Modern"),
            tickfont=(10, "Computer Modern"),
            guidefont=(14, "Computer Modern"),
            legend_font_pointsize = 14
        )
        Plots.scatter!(ps, x, u[:, i], label="",
        markersize=3)
    end
    Plots.plot!(
        titlefont=(14, "Computer Modern"),
        tickfont=(10, "Computer Modern"),
        guidefont=(14, "Computer Modern"),
        xlabel=L"x",
        ylabel=L"u_h(x, t)"
    )
    println("done!")

    println("Saving Samarski...")
    Plots.savefig(ps, joinpath("output/", "problem_2_implicit_Samarski.pdf"))


    return nothing
end


function high_precision_numerics()
    function analytical_solution(x, t)
        σ = 2
        ϰ₀ = 0.5
        x₁ = 0
        c = 5
        if x ≥ x₁ + c*(t + 0.1)
            return 0
        else
            return (σ*c*ϰ₀^(-1) * (c*(t + 0.1) + x₁ - x))^(1 / σ)
        end
    end

    k(u) = 0.5 * u^2
    f(x, t) = 0
    mu1(t) = 10 * √(t + 0.1)
    mu2(t) = 0
    u0(x) = (x < 0.5) ? (2 * √(5 * (0.5 - x))) : 0

    T = 0.1
    Nx = 5001
    Nt = 5001
    x = range(0, 1, length = Nx)
    t = range(0, T, length = Nt)
    number_of_t_dots = 50
    number_of_x_dots = 50
    stepik_t = Nt ÷ number_of_t_dots
    stepik_x = Nx ÷ number_of_x_dots

    println("Solving numerically...")
    u_h = implicit_scheme(k, f, u0, mu1, mu2, T, Nx, Nt)
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


    println("Plotting surface solution")
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
        camera=(-15, 20),
        colormap_name = "viridis",
        extra_kwargs =:subplot,
        cbar=false
    )
    println("done!")

    println("Saving surface...")
    Plots.savefig(p, joinpath("output/", "problem_2_implicit_surface_high.pdf"))
    println("done!")

    println("Plotting error surface...")
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
        colormap_name = "viridis",
        extra_kwargs =:subplot,
        camera=(15, 20),
        cbar=false
    )
    println("done!")

    println("Saving...")
    Plots.savefig(p, joinpath("output/", "problem_2_implicit_error_surface_high.pdf"))
    println("done!")

    return nothing
end


function error_Samarski()

    function analytical_solution(x, t)
        σ = 2
        ϰ₀ = 0.5
        x₁ = 0
        c = 5
        if x ≥ x₁ + c*(t + 0.1)
            return 0
        else
            return (σ*c*ϰ₀^(-1) * (c*(t + 0.1) + x₁ - x))^(1 / σ)
        end
    end

    k(u) = 0.5 * u^2
    f(x, t) = 0
    mu1(t) = 10 * √(t + 0.1)
    mu2(t) = 0
    u0(x) = (x < 0.5) ? (2 * √(5 * (0.5 - x))) : 0

    T = 0.1
    Nx = 51
    Nt = 501
    x = range(0, 1, length = Nx)
    t = range(0, T, length = Nt)

    println("Solving numerically...")
    u_h = implicit_scheme(k, f, u0, mu1, mu2, T, Nx, Nt)
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
    index = [1 101 301 501]
    xlimits = (0, 1)
    ps = Plots.plot()
    for i ∈ index
        Plots.plot!(
            ps, x, err[:, i], xlims=xlimits, label="",
            titlefont=(14, "Computer Modern"),
            tickfont=(10, "Computer Modern"),
            guidefont=(14, "Computer Modern"),
            legend_font_pointsize = 14
        )
        # Plots.scatter!(ps, x, u[:, i], label="",
        # markersize=3)
    end
    Plots.plot!(
        titlefont=(14, "Computer Modern"),
        tickfont=(10, "Computer Modern"),
        guidefont=(14, "Computer Modern"),
        xlabel=L"x",
        ylabel=L"|u(x, t) - u_h(x, t)|"
    )
    println("done!")

    println("Saving Samarski...")
    Plots.savefig(ps, joinpath("output/", "problem_2_implicit_error_Samarski.pdf"))

    return nothing
end


function error_surface_2()
    function analytical_solution(x, t)
        σ = 2
        ϰ₀ = 0.5
        x₁ = 0
        c = 5
        if x ≥ x₁ + c*(t + 0.1)
            return 0
        else
            return (σ*c*ϰ₀^(-1) * (c*(t + 0.1) + x₁ - x))^(1 / σ)
        end
    end

    k(u) = 0.5 * u^2
    f(x, t) = 0
    mu1(t) = 10 * √(t + 0.1)
    mu2(t) = 0
    u0(x) = (x < 0.5) ? (2 * √(5 * (0.5 - x))) : 0

    T = 0.1
    Nx = 51
    Nt = 501
    x = range(0, 1, length = Nx)
    t = range(0, T, length = Nt)

    println("Solving numerically...")
    u_h = implicit_scheme(k, f, u0, mu1, mu2, T, Nx, Nt)
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
    p = Plots.plot(
        x, t .+ 0.1, err',
        label="",
        titlefont=(14, "Computer Modern"),
        tickfont=(10, "Computer Modern"),
        guidefont=(14, "Computer Modern"),
        xlabel=L"x",
        ylabel=L"t",
        zlabel=L"|u(x, t) - u_h(x, t)|",
        st = :surface,
        colormap_name = "viridis",
        camera=(15, 20),
        cbar=false
    )
    println("done!")

    println("Saving...")
    Plots.savefig(p, joinpath("output/", "problem_2_implicit_error_surface.pdf"))
    println("done!")

    return nothing
end