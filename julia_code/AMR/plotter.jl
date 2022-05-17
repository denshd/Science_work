using Plots, LaTeXStrings, Printf

include("solver.jl")

function plot_gif(
    solution::UniformGridSolution;
    t1 = solution.time_grid.t[1],
    t2 = solution.time_grid.t[end],
    color = :black,
    number_of_frames = min(100, solution.time_grid.Ntime),
    fps = 60,
    filename = "solution.gif",
    picture_size=(1280, 1024)
    )

    time = solution.time_grid.t
    Ntime = solution.time_grid.Ntime
    x = solution.spacial_grid.cells
    u = solution.u

    time_step = div(Ntime, number_of_frames)

    xlimits = (x[1], x[end])
    ylimits = (min(u...), max(u...))

    anim = Animation()
    for i = 1:time_step:Ntime
        label = latexstring(@sprintf "t = %.3f" time[i])
        p = plot(
            x, u[:, i],
            color=color,
            xlims=xlimits,
            ylims=ylimits,
            title=L"solution $u(x, t)$",
            xlabel=L"x", ylabel=L"u",
            label=label,
            size=picture_size
        )
        frame(anim, p)
    end

    G = gif(anim, filename);
end