include("schemes.jl")

using Plots, GLMakie, LaTeXStrings, LinearAlgebra, BenchmarkTools, Printf
pgfplotsx()


function analytical_solution(x, y, t)
    if (t > (x + 2y))
        return 0.5 * √(-1 + √(1 + 16 * (t - x - 2y)))
    else
        return 0
    end
end


function find_max(X::Array{Cdouble, 3})
    temp_max = -Inf
    for x in X
        if (x > temp_max)
            temp_max = x
        end
    end
    return temp_max
end


function find_min(X::Array{Cdouble, 3})
    temp_min = +Inf
    for x in X
        if (x < temp_min)
            temp_min = x
        end
    end
    return temp_min
end


function main()
    println("Prepearing...")
    f(x, y, t) = 0.0
    k1(u) = 4 * u^4
    k2(u) = 0.25 * u^2
    u0(x, y) = 0.0

    mu_1(y, t) = (t > 2y) ? (0.5√(-1 + √(1 + 16(t - 2y)))) : 0
    mu1(y, t) = (t > 30 + 2y) ? (0.5√(-1 + √(1 + 16(t - 30 - 2y)))) : 0
    mu_2(x, t) = (t > x) ? (0.5√(-1 + √(1 + 16(t - x)))) : 0
    mu2(x, t) = (t > x + 40) ? (0.5√(-1 + √(1 + 16(t - x - 40)))) : 0

    Lx = 30
    Ly = 20
    T = 50

    Nx = 501
    Ny = 501
    Nt = 1001

    x = range(0, Lx, length = Nx)
    dx = step(x)
    y = range(0, Ly, length = Ny)
    dy = step(y)
    t = range(0, T, length = Nt)
    dt = step(t)

    println("done!")

    println("Solving...")
    local_scheme(k1, k2, f, u0, mu_1, mu1, mu_2, mu2, 10, 10, 50, 10, 10, 10);
    println("compiled, now real...")
    @time u = local_scheme(k1, k2, f, u0, mu_1, mu1, mu_2, mu2, Lx, Ly, T, Nx, Ny, Nt)
    println("done!")

    println("Plotting...")
    number_of_t_dots = 20
    number_of_x_dots = 20
    number_of_y_dots = 20
    stepik_t = Nt ÷ min(number_of_t_dots, Nt)
    stepik_x = Nx ÷ min(number_of_x_dots, Nx)
    stepik_y = Ny ÷ min(number_of_y_dots, Ny)
    number_of_plots = 9
    time_step = Nt ÷ min(number_of_plots, Nt)

    limits = (
        (x[1], x[end]),
        (y[1], y[end]),
        (find_min(u), find_max(u))
    )

    println("Creating news...")
    new_x = collect(range(0, Lx, number_of_x_dots))
    sx = step(range(0, Lx, number_of_x_dots))

    new_y = collect(range(0, Ly, number_of_y_dots))
    sy =step(range(0, Ly, number_of_y_dots))

    new_t = collect(range(0, T, number_of_plots))
    st = step(range(0, T, number_of_plots))

    new_u = zeros(length(new_x), length(new_y), length(new_t))
    for i in 1:length(new_x)
        for k in 1:length(new_y)
            for j in 1:length(new_t)
                new_u[i, k, j] = u[floor(Int, 1 + (i - 1)*sx), floor(Int, 1 + (k - 1)*sy), floor(Int, 1 + (j - 1)*st)]
            end
        end
    end
    println("done!")

    p = []

    for i in 1:length(new_t)
        push!(p, Plots.plot(
            new_x, new_y, new_u[:, :, i]',
            st = :wireframe,
            label="",
            titlefont=(14, "Computer Modern"),
            tickfont=(10, "Computer Modern"),
            guidefont=(14, "Computer Modern"),
            xlabel=L"x",
            ylabel=L"y",
            zlabel=L"u_h(x, t)",
            cbar=false,
            camera=(120, 20),
            title = latexstring(@sprintf("t = %.3f", new_t[i]))
        ))

    # for i in 1:number_of_plots
    #     push!(p, Plots.plot(
    #         x[1:stepik_x:end], y[1:stepik_y:end], u[1:stepik_x:end, 1:stepik_y:end, (i - 1) * time_step + 1]',
    #         st = :wireframe,
    #         label="",
    #         titlefont=(14, "Computer Modern"),
    #         tickfont=(10, "Computer Modern"),
    #         guidefont=(14, "Computer Modern"),
    #         xlabel=L"x",
    #         ylabel=L"y",
    #         zlabel=L"u_h(x, t)",
    #         cbar=false,
    #         camera=(120, 30),
    #         title = latexstring(@sprintf("t = %.3f", t[(i - 1)*time_step + 1]))
    #     ))
        # push!(p, Plots.heatmap(
        #     x, y, u[:, :, (i - 1)*time_step + 1],
        #     label="",
        #     titlefont=(14, "Computer Modern"),
        #     tickfont=(10, "Computer Modern"),
        #     guidefont=(14, "Computer Modern"),
        #     cbar_lims = 
        #     title = latexstring(@sprintf("%.3f", t[(i - 1)*time_step + 1]))
        # ))
    end
    res = Plots.plot(p..., size=(800, 1000))
    println("done!")

    println("Saving...")
    Plots.savefig(res, joinpath("output/", "problem_3_loc_wireframe.pdf"))
    println("done!")
end