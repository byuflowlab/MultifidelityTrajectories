function plot_states(ts, xs, nx; figure_name = "states", clear_figure = true, save_figure = true, fig_name = "states.png", color = "blue", stl = "-", label = "")
    fig_states = PP.figure(figure_name)
    if clear_figure; fig_states.clear(); end

    fig_states.add_subplot(231, ylabel = LS.L"x[m]")
    fig_states.add_subplot(234, ylabel = LS.L"y[m]", xlabel = LS.L"t[s]")
    fig_states.add_subplot(232, ylabel = LS.L"v_x[m/s]")
    fig_states.add_subplot(235, ylabel = LS.L"v_y[m/s]", xlabel = LS.L"t[s]")
    fig_states.add_subplot(233, ylabel = LS.L"\theta[^\circ]")
    fig_states.add_subplot(236, ylabel = LS.L"\dot{\theta}[^\circ/s]", xlabel = LS.L"t[s]")

    axs = fig_states.get_axes()
    for (i,ax) in enumerate(axs)
        ax.plot(ts, index_states(xs, nx, i), color = color, stl, label = label)
    end

    fig_states.tight_layout()
    fig_states.set_size_inches(11, 7, forward=true)

    if save_figure
        if !isdir(joinpath(topdirectory, "data", "plots", TODAY)); mkpath(joinpath(topdirectory, "data", "plots", TODAY)); end
        fig_states.savefig(joinpath(topdirectory, "data", "plots", TODAY, fig_name))
    end
end
