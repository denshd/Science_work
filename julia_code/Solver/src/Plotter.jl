using Plots, LaTeXStrings

export plot_gif, plot_grid

function plot_gif(sol::Solution{<:Grid1};
    t1 = sol.g.t[1],
    t2 = sol.g.t[2],
    color="black",
    number_of_frames = min(100, length(sol.g.t)),
    fps = 30,
    filename = "solution.gif",
    picture_size=(1280, 1024)
    )

    time_step = div(length(sol.g.t), number_of_frames)

    xlimits = (sol.g.x[1], sol.g.x[end])
    ylimits = (min(sol.u...), max(sol.u...))

    anim = Animation()
    for i = 1 : time_step : sol.g.Nt
        p = plot(
            sol.g.x, sol.u[:, i],
            color=color,
            xlims=xlimits,
            ylims=ylimits,
            title=L"solution $u(x, t)$",
            xlabel=L"x", ylabel=L"u",
            label="t=$(round(sol.g.t[i], digits=3))",
            size=picture_size
        )
        frame(anim, p)
    end

    G = gif(anim, filename)
end


function plot_gif(sol::Solution2{<:Grid2};
    t1 = sol.g.t[1],
    t2 = sol.g.t[2],
    number_of_frames = min(100, length(sol.g.t)),
    fps = 30,
    filename = "heatmap.gif",
    picture_size=(1280, 1024)
    )

    time_step = div(length(sol.g.t), number_of_frames)

    limits = (
        (sol.g.x[1][1], sol.g.x[1][end]),
        (sol.g.x[2][1], sol.g.x[2][end]),
        (min(sol.u...), max(sol.u...))
    )

    anim = Animation()
    for i = 1 : time_step : sol.g.Nt
        heatmap(sol.g.x[1], sol.g.x[2], sol.u[:, :, i], c = :thermal, size=picture_size)
        frame(anim)
    end

    G = gif(anim, filename)
end


function plot_grid(g::Grid1; filename="grid.pdf")
    p = plot(scatter([(g.x[i], g.t[j]) for i in 1:g.Nx for j in 1:g.Nt]))
    savefig(p, filename)
end


function plot_grid(g::Grid2; filename="grid.pdf")
    p = plot(scatter([(g.x[1][i], g.x[2][j]) for i in 1:g.N[1] for j in 1:g.N[2]]))
    savefig(p, filename)
end