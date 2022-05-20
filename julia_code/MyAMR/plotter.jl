using Plots, LaTeXStrings, Printf

include("grids.jl")

"""
TODO: Всё это изменить под новые структуры, написать визуализатор Level'ов
"""


function plot_grid(
    level::Level{T};
    filename=joinpath("output", "grid_test.pdf")
    ) where T

    p = plot()
    plot_grid!(p, level)

    savefig(p, filename)
end


function plot_grid!(
    p,
    level::Level{T}
    ) where T

    for block in level.blocks
        plot_grid!(p, block, level.level_number)
    end

    if has_sublevel(level)
        plot_grid!(p, get_sublevel(level))
    end
end


function plot_grid!(
    p,
    block::Block{T},
    level_number::Integer
    ) where T

    plot_grid!(p, block.spacial_grid, level_number)
end


function plot_grid!(
    p,
    g::UniformGrid2{T},
    level_number::Integer
    ) where T

    X = g.g[1].x
    Y = g.g[2].x

    for x in X
        plot!(
            p,
            [x, x], [Y[1], Y[end]],
            label="",
            color=:black,
            linewidth=(1 / (2)^level_number)
        )
    end
    for y in Y
        plot!(
            p,
            [X[1], X[end]], [y, y],
            label="",
            color=:black,
            linewidth=(1 / (2)^level_number)
        )
    end
end


function plot_gif(
    solution::Solution;
    t1 = solution.time_grid.x[1],
    t2 = solution.time_grid.x[end],
    number_of_frames = min(100, solution.time_grid.N),
    fps = 60,
    filename = "solution.gif",
    picture_size=(1280, 1024)
    )

    time = solution.time_grid.x
    Ntime = solution.time_grid.N
    x = solution.levels[1].sublevel.sublevel.blocks[6].spacial_grid.g[1].x
    y = solution.levels[1].sublevel.sublevel.blocks[6].spacial_grid.g[2].x
    # cbar_limits = (min(u...), max(u...))

    time_step = div(Ntime, number_of_frames)

    # limits = (
    #     (x[1], x[end]),
    #     (y[1], y[end]),
    #     (min(u...), max(u...))
    # )

    anim = Animation()
    for i = 1:time_step:Ntime
        label = latexstring(@sprintf "t = %.3f" time[i])
        heatmap(
            x, y,
            solution.levels[i].sublevel.sublevel.blocks[6].u_new',
            c = :thermal,
            # cbar_lims=cbar_limits,
            size = picture_size
        )
        frame(anim)
    end

    G = gif(anim, filename, fps=1);
end

# """
# Строит анимацию решение одномерного уравнения теплопроводности
# """
# function plot_gif(
#     solution::UniformGridSolution;
#     t1 = solution.time_grid.t[1],
#     t2 = solution.time_grid.t[end],
#     color = :black,
#     number_of_frames = min(100, solution.time_grid.Ntime),
#     fps = 60,
#     filename = "solution.gif",
#     picture_size=(1280, 1024)
#     )

#     time = solution.time_grid.t
#     Ntime = solution.time_grid.Ntime
#     x = solution.spacial_grid.cells
#     u = solution.u

#     time_step = div(Ntime, number_of_frames)

#     xlimits = (x[1], x[end])
#     ylimits = (min(u...), max(u...))

#     anim = Animation()
#     for i = 1:time_step:Ntime
#         label = latexstring(@sprintf "t = %.3f" time[i])
#         p = plot(
#             x, u[:, i],
#             color=color,
#             xlims=xlimits,
#             ylims=ylimits,
#             title=L"solution $u(x, t)$",
#             xlabel=L"x", ylabel=L"u",
#             label=label,
#             size=picture_size
#         )
#         frame(anim, p)
#     end

#     G = gif(anim, filename);
# end


# """
# Строит анимацию решения двумерного уравнения теплопроводности
# """
# function plot_gif(
#     solution::UniformGrid2Solution;
#     t1 = solution.time_grid.t[1],
#     t2 = solution.time_grid.t[end],
#     color = :black,
#     number_of_frames = min(100, solution.time_grid.Ntime),
#     fps = 60,
#     filename = "solution.gif",
#     picture_size=(1280, 1024)
#     )

#     time = solution.time_grid.t
#     Ntime = solution.time_grid.Ntime
#     x = solution.spacial_grid.g[1].cells
#     y = solution.spacial_grid.g[2].cells
#     u = solution.u
#     cbar_limits = (min(u...), max(u...))

#     time_step = div(Ntime, number_of_frames)

#     limits = (
#         (x[1], x[end]),
#         (y[1], y[end]),
#         (min(u...), max(u...))
#     )

#     anim = Animation()
#     for i = 1:time_step:Ntime
#         label = latexstring(@sprintf "t = %.3f" time[i])
#         heatmap(
#             x, y, u[:, :, i]',
#             c = :thermal,
#             cbar_lims=cbar_limits,
#             size = picture_size
#         )
#         frame(anim)
#     end

#     G = gif(anim, filename);
# end


# """
# Строит одномерную сетку
# """
# function plot_grid(g::UniformGrid; filename="grid.pdf")

#     # Построение точечек на пересечении линий
#     p = scatter(
#         [(g.faces[i], g.faces[j]) for i in 1:g.Nfaces for j in 1:2],
#         label="",
#         color=:black,
#         markersize=3,
#         markeralpha=0.5,
#         markerstrokecolor=:black
#     )

#     # Построение линий (faces)
#     for i = 1:g.Nfaces
#         plot!(
#             [g.faces[i], g.faces[i]], [g.faces[1], g.faces[2]],
#             label="",
#             color=:black
#         )
#     end
#     plot!(
#         [g.faces[1], g.faces[end]], [g.faces[1], g.faces[1]],
#         label="",
#         color=:black
#     )
#     plot!(
#         [g.faces[1], g.faces[end]], [g.faces[2], g.faces[2]],
#         label="",
#         color=:black
#     )

#     # Построение кружочков
#     scatter!(
#         [(g.cells[i], g.cells[1]) for i in 1:g.Ncells],
#         label="",
#         color=:white,
#         markersize=3,
#         markerstrokecolor=:red
#     )
#     savefig(p, filename);
# end


# """
# Строит двумерную сетку
# """
# function plot_grid(g::UniformGrid2; filename="grid.pdf")

#     # Построение точечек
#     p = scatter(
#         [(g.g[1].faces[i], g.g[2].faces[j]) for i in 1:g.g[1].Nfaces for j in 1:g.g[2].Nfaces],
#         label="",
#         color=:black,
#         markersize=3,
#         markerstrokecolor=:black
#     )

#     # Построение линий
#     for i in 1:g.g[1].Nfaces
#         plot!(
#             [g.g[1].faces[i], g.g[1].faces[i]], [g.g[2].faces[1], g.g[2].faces[end]],
#             label="",
#             color=:black
#         )
#     end
#     for i in 1:g.g[2].Nfaces
#         plot!(
#             [g.g[1].faces[1], g.g[1].faces[end]], [g.g[2].faces[i], g.g[2].faces[i]],
#             label="",
#             color=:black
#         )
#     end

#     # Построение кружочков внутри клеточек
#     scatter!(
#         [(x, y) for x in g.g[1].cells for y in g.g[2].cells],
#         label="",
#         color=:white,
#         markersize=3,
#         markerstrokecolor=:red
#     )

#     savefig(p, filename);
# end