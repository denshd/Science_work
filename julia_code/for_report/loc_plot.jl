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
    println("Preparing...")
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

    Nx = 31
    Ny = 31
    Nt = 101

    x = range(0, Lx, length = Nx); dx = step(x)
    y = range(0, Ly, length = Ny); dy = step(y)
    t = range(0, T, length = Nt); dt = step(t)
    println("done!")

    println("Solving...")
    local_scheme(k1, k2, f, u0, mu_1, mu1, mu_2, mu2, 10, 10, 50, 10, 10, 10);
    @time u_h = local_scheme(k1, k2, f, u0, mu_1, mu1, mu_2, mu2, Lx, Ly, T, Nx, Ny, Nt)
    println("done!")

    println("Calculating analytic solution and error...")
    u = zeros(size(u_h))
    err = zeros(size(u))
    for i in 1:Nx
        for k in 1:Ny
            for j in 1:Nt
                u[i, k, j] = analytical_solution(x[i], y[k], t[j])
                err[i, k, j] = abs(u[i, k, j] - u_h[i, k, j])
            end
        end
    end
    println("done!")

    println("Prepare for plotting...")
    limits = ((x[1], x[end]), (y[1], y[end]), (find_min(u_h), find_max(u_h)), (find_min(err), find_max(err)))
    number_of_plots = 9
    time_step = Nt ÷ min(Nt, number_of_plots)
    println("done!")

    println("Plotting...")
    p_h = []
    p_err = []

    for i in 1:number_of_plots
        push!(p_h, Plots.plot(
                x, y, u_h[:, :, (i - 1) * time_step + 1]',
                st = :wireframe,
                label = "",
                titlefont=(14, "Computer Modern"),
                tickfont=(10, "Computer Modern"),
                guidefont=(14, "Computer Modern"),
                xlabel=L"x",
                ylabel=L"y",
                zlabel="",
                cbar=false,
                camera=(120, 20),
                title = latexstring(@sprintf("t = %.2f", t[(i - 1) * time_step + 1])),
                xlims=limits[1],
                ylims=limits[2],
                zlims=limits[3]
            )
        )
        push!(p_err, Plots.plot(
                x, y, err[:, :, (i - 1) * time_step + 1]',
                st = :wireframe,
                label = "",
                titlefont=(14, "Computer Modern"),
                tickfont=(10, "Computer Modern"),
                guidefont=(14, "Computer Modern"),
                xlabel=L"x",
                ylabel=L"y",
                zlabel="",
                cbar=false,
                camera=(70, 20),
                title = latexstring(@sprintf("t = %.3f", t[(i - 1) * time_step + 1])),
                xlims=limits[1],
                ylims=limits[2],
                zlims=limits[4]
            )
        )
    end

    res_h = Plots.plot(p_h..., size=(800, 1000))
    res_err = Plots.plot(p_err..., size=(800, 1000))
    println("done!")

    println("Saving...")
    Plots.savefig(res_h, joinpath("output/", "problem_3_loc_wireframe.pdf"))
    Plots.savefig(res_err, joinpath("output/", "problem_3_loc_err_wireframe.pdf"))
    println("done!")
end



# function temp_plot()
#     println("Prepearing...")
#     f(x, y, t) = 0.0
#     k1(u) = 4 * u^4
#     k2(u) = 0.25 * u^2
#     u0(x, y) = 0.0

#     mu_1(y, t) = (t > 2y) ? (0.5√(-1 + √(1 + 16(t - 2y)))) : 0
#     mu1(y, t) = (t > 30 + 2y) ? (0.5√(-1 + √(1 + 16(t - 30 - 2y)))) : 0
#     mu_2(x, t) = (t > x) ? (0.5√(-1 + √(1 + 16(t - x)))) : 0
#     mu2(x, t) = (t > x + 40) ? (0.5√(-1 + √(1 + 16(t - x - 40)))) : 0

#     Lx = 30
#     Ly = 20
#     T = 50

#     Nx = 21
#     Ny = 21
#     Nt = 101

#     x = range(0, Lx, length = Nx)
#     dx = step(x)
#     y = range(0, Ly, length = Ny)
#     dy = step(y)
#     t = range(0, T, length = Nt)
#     dt = step(t)

#     println("done!")

#     println("Solving...")
#     local_scheme(k1, k2, f, u0, mu_1, mu1, mu_2, mu2, 10, 10, 50, 10, 10, 10);
#     println("compiled, now real...")
#     @time u_h = local_scheme(k1, k2, f, u0, mu_1, mu1, mu_2, mu2, Lx, Ly, T, Nx, Ny, Nt)
#     println("done!")

#     println("Plotting...")


#     println("Calculating analytic solution and error...")
#     u = zeros(size(u_h))
#     err = zeros(size(u))
#     for i in 1:length(x)
#         for k in 1:length(y)
#             for j in 1:length(t)
#                 u[i, k, j] = analytical_solution(x[i], y[k], t[j])
#                 err[i, k, j] = abs(u[i, k, j] - u_h[i, k, j])
#                 # err[i, k, j] = u_h[i, k, j]
#             end
#         end
#     end

#     limits = (
#         (x[1], x[end]),
#         (y[1], y[end]),
#         (find_min(err), find_max(err))
#     )

#     p = []

#     number_of_plots = 9
#     time_step = Nt ÷ min(number_of_plots, Nt)
#     for i in 1:number_of_plots
#         push!(p, Plots.plot(
#                 x, y, err[:, :, (i - 1) * time_step + 1]',
#                 st = :wireframe,
#                 label = "",
#                 titlefont=(14, "Computer Modern"),
#                 tickfont=(10, "Computer Modern"),
#                 guidefont=(14, "Computer Modern"),
#                 xlabel=L"x",
#                 ylabel=L"y",
#                 zlabel=L"|u_h(x, t) - u(x, t)|",
#                 cbar=false,
#                 camera=(70, 20),
#                 title = latexstring(@sprintf("t = %.3f", t[(i - 1) * time_step + 1]))
#             )
#         )
#     end


#     res = Plots.plot(p..., size=(800, 1000))
#     println("done!")

#     println("Saving...")
#     Plots.savefig(res, joinpath("output/", "problem_3_loc_err_wireframe.pdf"))
#     println("done!")
# end


# function error_plot()
#     println("Prepearing...")
#     f(x, y, t) = 0.0
#     k1(u) = 4 * u^4
#     k2(u) = 0.25 * u^2
#     u0(x, y) = 0.0

#     mu_1(y, t) = (t > 2y) ? (0.5√(-1 + √(1 + 16(t - 2y)))) : 0
#     mu1(y, t) = (t > 30 + 2y) ? (0.5√(-1 + √(1 + 16(t - 30 - 2y)))) : 0
#     mu_2(x, t) = (t > x) ? (0.5√(-1 + √(1 + 16(t - x)))) : 0
#     mu2(x, t) = (t > x + 40) ? (0.5√(-1 + √(1 + 16(t - x - 40)))) : 0

#     Lx = 30
#     Ly = 20
#     T = 50

#     Nx = 21
#     Ny = 21
#     Nt = 101

#     x = range(0, Lx, length = Nx)
#     dx = step(x)
#     y = range(0, Ly, length = Ny)
#     dy = step(y)
#     t = range(0, T, length = Nt)
#     dt = step(t)

#     println("done!")

#     println("Solving...")
#     local_scheme(k1, k2, f, u0, mu_1, mu1, mu_2, mu2, 10, 10, 50, 10, 10, 10);
#     println("compiled, now real...")
#     @time u_h = local_scheme(k1, k2, f, u0, mu_1, mu1, mu_2, mu2, Lx, Ly, T, Nx, Ny, Nt)
#     println("done!")

#     println("Plotting...")
#     number_of_t_dots = 51
#     number_of_x_dots = 51
#     number_of_y_dots = 21
#     stepik_t = Nt ÷ min(number_of_t_dots, Nt)
#     stepik_x = Nx ÷ min(number_of_x_dots, Nx)
#     stepik_y = Ny ÷ min(number_of_y_dots, Ny)
#     number_of_plots = 9
#     time_step = Nt ÷ min(number_of_plots, Nt)

#     limits = (
#         (x[1], x[end]),
#         (y[1], y[end]),
#         (find_min(u_h), find_max(u_h))
#     )

#     println("Creating news...")
#     new_x = collect(range(0, Lx, number_of_x_dots))
#     sx = step(range(0, Lx, number_of_x_dots))

#     new_y = collect(range(0, Ly, number_of_y_dots))
#     sy =step(range(0, Ly, number_of_y_dots))

#     new_t = collect(range(0, T, number_of_plots))
#     st = step(range(0, T, number_of_plots))

#     new_u = zeros(length(new_x), length(new_y), length(new_t))
#     for i in 1:length(new_x)
#         for k in 1:length(new_y)
#             for j in 1:length(new_t)
#                 new_u[i, k, j] = u_h[floor(Int, 1 + (i - 1)*sx), floor(Int, 1 + (k - 1)*sy), floor(Int, 1 + (j - 1)*st)]
#             end
#         end
#     end
#     println("done!")


#     println("Calculating analytic solution and error...")
#     u = zeros(size(new_u))
#     err = zeros(size(u))
#     for i in 1:length(new_x)
#         for k in 1:length(new_y)
#             for j in 1:length(new_t)
#                 u[i, k, j] = analytical_solution(new_x[i], new_y[k], new_t[j])
#                 err[i, k, j] = abs(u[i, k, j] - new_u[i, k, j])
#                 # err[i, k, j] = u[i, k, j]
#             end
#         end
#     end

#     p = []

#     for i in 1:length(new_t)
#         push!(p, Plots.plot(
#             new_x, new_y, err[:, :, i]',
#             st = :wireframe,
#             label="",
#             titlefont=(14, "Computer Modern"),
#             tickfont=(10, "Computer Modern"),
#             guidefont=(14, "Computer Modern"),
#             xlabel=L"x",
#             ylabel=L"y",
#             zlabel=L"u_h(x, t)",
#             cbar=false,
#             camera=(120, 20),
#             title = latexstring(@sprintf("t = %.3f", new_t[i]))
#         ))

#     # for i in 1:number_of_plots
#     #     push!(p, Plots.plot(
#     #         x[1:stepik_x:end], y[1:stepik_y:end], u[1:stepik_x:end, 1:stepik_y:end, (i - 1) * time_step + 1]',
#     #         st = :wireframe,
#     #         label="",
#     #         titlefont=(14, "Computer Modern"),
#     #         tickfont=(10, "Computer Modern"),
#     #         guidefont=(14, "Computer Modern"),
#     #         xlabel=L"x",
#     #         ylabel=L"y",
#     #         zlabel=L"u_h(x, t)",
#     #         cbar=false,
#     #         camera=(120, 30),
#     #         title = latexstring(@sprintf("t = %.3f", t[(i - 1)*time_step + 1]))
#     #     ))
#         # push!(p, Plots.heatmap(
#         #     x, y, u[:, :, (i - 1)*time_step + 1],
#         #     label="",
#         #     titlefont=(14, "Computer Modern"),
#         #     tickfont=(10, "Computer Modern"),
#         #     guidefont=(14, "Computer Modern"),
#         #     cbar_lims = 
#         #     title = latexstring(@sprintf("%.3f", t[(i - 1)*time_step + 1]))
#         # ))
#     end
#     res = Plots.plot(p..., size=(800, 1000))
#     println("done!")

#     println("Saving...")
#     Plots.savefig(res, joinpath("output/", "problem_3_loc_err_wireframe.pdf"))
#     println("done!")
# end


# function main()
#     println("Prepearing...")
#     f(x, y, t) = 0.0
#     k1(u) = 4 * u^4
#     k2(u) = 0.25 * u^2
#     u0(x, y) = 0.0

#     mu_1(y, t) = (t > 2y) ? (0.5√(-1 + √(1 + 16(t - 2y)))) : 0
#     mu1(y, t) = (t > 30 + 2y) ? (0.5√(-1 + √(1 + 16(t - 30 - 2y)))) : 0
#     mu_2(x, t) = (t > x) ? (0.5√(-1 + √(1 + 16(t - x)))) : 0
#     mu2(x, t) = (t > x + 40) ? (0.5√(-1 + √(1 + 16(t - x - 40)))) : 0

#     Lx = 30
#     Ly = 20
#     T = 50

#     Nx = 501
#     Ny = 501
#     Nt = 1001

#     x = range(0, Lx, length = Nx)
#     dx = step(x)
#     y = range(0, Ly, length = Ny)
#     dy = step(y)
#     t = range(0, T, length = Nt)
#     dt = step(t)

#     println("done!")

#     println("Solving...")
#     local_scheme(k1, k2, f, u0, mu_1, mu1, mu_2, mu2, 10, 10, 50, 10, 10, 10);
#     println("compiled, now real...")
#     @time u = local_scheme(k1, k2, f, u0, mu_1, mu1, mu_2, mu2, Lx, Ly, T, Nx, Ny, Nt)
#     println("done!")

#     println("Plotting...")
#     number_of_x_dots = 21
#     number_of_y_dots = 21
#     number_of_plots = 9

#     println("Creating news...")
#     new_x = collect(range(0, Lx, number_of_x_dots))
#     sx = step(range(0, Lx, number_of_x_dots))

#     new_y = collect(range(0, Ly, number_of_y_dots))
#     sy =step(range(0, Ly, number_of_y_dots))

#     new_t = collect(range(0, T, number_of_plots))
#     st = step(range(0, T, number_of_plots))

#     new_u = zeros(length(new_x), length(new_y), length(new_t))
#     for i in 1:length(new_x)
#         for k in 1:length(new_y)
#             for j in 1:length(new_t)
#                 new_u[i, k, j] = u[floor(Int, 1 + (i - 1)*sx), floor(Int, 1 + (k - 1)*sy), floor(Int, 1 + (j - 1)*st)]
#             end
#         end
#     end

#     limits = (
#         (new_x[1], new_x[end]),
#         (new_y[1], new_y[end]),
#         (find_min(u), find_max(u))
#     )
#     println("done!")

#     p = []

#     for i in 1:length(new_t)
#         push!(p, Plots.plot(
#             new_x, new_y, new_u[:, :, i]',
#             st = :wireframe,
#             label="",
#             titlefont=(14, "Computer Modern"),
#             tickfont=(10, "Computer Modern"),
#             guidefont=(14, "Computer Modern"),
#             xlabel=L"x",
#             ylabel=L"y",
#             zlabel=L"u_h(x, t)",
#             xlims=limits[1],
#             ylims=limits[2],
#             zlims=limits[3],
#             cbar=false,
#             camera=(120, 20),
#             title = latexstring(@sprintf("t = %.3f", new_t[i]))
#         ))

#     # for i in 1:number_of_plots
#     #     push!(p, Plots.plot(
#     #         x[1:stepik_x:end], y[1:stepik_y:end], u[1:stepik_x:end, 1:stepik_y:end, (i - 1) * time_step + 1]',
#     #         st = :wireframe,
#     #         label="",
#     #         titlefont=(14, "Computer Modern"),
#     #         tickfont=(10, "Computer Modern"),
#     #         guidefont=(14, "Computer Modern"),
#     #         xlabel=L"x",
#     #         ylabel=L"y",
#     #         zlabel=L"u_h(x, t)",
#     #         cbar=false,
#     #         camera=(120, 30),
#     #         title = latexstring(@sprintf("t = %.3f", t[(i - 1)*time_step + 1]))
#     #     ))
#         # push!(p, Plots.heatmap(
#         #     x, y, u[:, :, (i - 1)*time_step + 1],
#         #     label="",
#         #     titlefont=(14, "Computer Modern"),
#         #     tickfont=(10, "Computer Modern"),
#         #     guidefont=(14, "Computer Modern"),
#         #     cbar_lims = 
#         #     title = latexstring(@sprintf("%.3f", t[(i - 1)*time_step + 1]))
#         # ))
#     end
#     res = Plots.plot(p..., size=(800, 1000))
#     println("done!")

#     println("Saving...")
#     Plots.savefig(res, joinpath("output/", "problem_3_loc_wireframe.pdf"))
#     println("done!")
# end