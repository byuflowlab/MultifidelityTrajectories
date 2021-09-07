#####
##### Model
#####
crc3, parameters, environment, alphas, stepsymbol, cl_stall = MT.build_crc3(; nwings)
crc3_model = MT.AircraftDynamics(crc3, parameters, environment, stepsymbol, cl_stall)

function X_2_Xdot(model, alpha, gamma; vinf = 1.0, rotors_on = false, U = zeros(length(crc3_model.aircraft.rotor_system.index)) * 2*pi/60.0)
    X = MT.steady_climb_state(alpha, gamma, vinf)
    Xdot = MT.dynamics(model, X, U; rotors_on)
    return X, Xdot
end

function test_xdot_cl_alpha(crc3_model;
    alphas = range(-40.0, stop=55.0, length=26) .* pi/180,
    gamma = 0 * pi/180,
    vinf = 1.0,
    rotors_on = false,
    U = zeros(length(crc3_model.aircraft.rotor_system.index)) * 2*pi/60.0,
    save_figures = true, fig_name = "test_dynamics.png", clear_figure = true, label_tag = ""
)
    # initialize vars
    Xs = zeros(7,length(alphas))
    Xdots = zeros(7,length(alphas))

    # build lift curves
    for (i,alpha) in enumerate(alphas)
        X, Xdot = X_2_Xdot(crc3_model, alpha, gamma; vinf, rotors_on, U)
        Xs[:,i] .= X
        Xdots[:,i] .= Xdot
    end

    # get forces to check (in global frame)
    _, CFs, CMs = test_cl_alpha(crc3_model; alphas, vinf, rotors_on, U, frame = MT.Global(), save_figures = false)

    # dimensionalize forces
    ρ = crc3_model.environment.ρ
    g = crc3_model.environment.g
    qinf = 0.5 * ρ * vinf^2
    m = crc3_model.aircraft.inertia_system.mass
    I = crc3_model.aircraft.inertia_system.inertia_y

    weight = MT.ST.@SVector [0.0, -g * m] # in global frame
    Fs = CFs .* (crc3_model.aircraft.wing_system.system.reference[1].S * qinf) .+ weight
    Ms = CMs .* (crc3_model.aircraft.wing_system.system.reference[1].S * crc3_model.aircraft.wing_system.system.reference[1].c * qinf)
    # predict dynamics: F = ma, T = Iα
    p_dots = deepcopy(Xs[3:4,:]) # in the body frame
    for (i,alpha) in enumerate(alphas)
        local theta = Xs[5,i]
        local R = MT.R_from_body(MT.Global(), alpha, theta)
        p_dots[:,i] .= R * p_dots[:,i] # rotate to global frame
    end
    v_dots = Fs ./ m # in the global frame
    for (i,alpha) in enumerate(alphas)
        local theta = Xs[5,i]
        local theta_dot = Xs[6,i]
        local R = MT.R_from_global(MT.Body(), alpha, theta)
        v_dots[:,i] .= R * v_dots[:,i] # rotate to body frame
            - theta_dot .* (MT.ST.@SArray [0. -1.; 1. 0.]) * Xs[3:4,i] # include coriolis force
    end
    θ_dots = Xs[6,:]
    θ_dot_dots = Ms ./ I
    Xdots_check = vcat(p_dots, v_dots, θ_dots', θ_dot_dots')

    fig_X_Xdot = PP.figure("test_Xdot")
    if clear_figure; fig_X_Xdot.clear(); end

    fig_X_Xdot.add_subplot(2, 3, 1, ylabel = MT.LS.L"P_x[m]")
    fig_X_Xdot.add_subplot(2, 3, 2, ylabel = MT.LS.L"V_x[m/s]")
    fig_X_Xdot.add_subplot(2, 3, 3, ylabel = MT.LS.L"\theta[^\circ]")
    fig_X_Xdot.add_subplot(2, 3, 4, ylabel = MT.LS.L"P_y[m]", xlabel = MT.LS.L"\alpha[^\circ]")
    fig_X_Xdot.add_subplot(2, 3, 5, ylabel = MT.LS.L"V_y[m/s]", xlabel = MT.LS.L"\alpha[^\circ]")
    fig_X_Xdot.add_subplot(2, 3, 6, ylabel = MT.LS.L"\dot{\theta}[^\circ/s]", xlabel = MT.LS.L"\alpha[^\circ]")

    axs_X_Xdot = fig_X_Xdot.get_axes()
    axs_X_Xdot[1].plot(alphas .* 180/pi, Xs[1,:], label = MT.LS.L"X" * label_tag)
    axs_X_Xdot[1].plot(alphas .* 180/pi, Xdots[1,:], label = MT.LS.L"\dot{X}" * label_tag)
    axs_X_Xdot[1].plot(alphas .* 180/pi, Xdots_check[1,:], label = "verify " * MT.LS.L"\dot{X}" * label_tag, "--")
    axs_X_Xdot[2].plot(alphas .* 180/pi, Xs[3,:], label = nothing)
    axs_X_Xdot[2].plot(alphas .* 180/pi, Xdots[3,:], label = nothing)
    axs_X_Xdot[2].plot(alphas .* 180/pi, Xdots_check[3,:], label = nothing, "--")
    axs_X_Xdot[3].plot(alphas .* 180/pi, Xs[5,:] .* 180/pi, label = nothing)
    axs_X_Xdot[3].plot(alphas .* 180/pi, Xdots[5,:] .* 180/pi, label = nothing)
    axs_X_Xdot[3].plot(alphas .* 180/pi, Xdots_check[5,:] .* 180/pi, label = nothing, "--")
    axs_X_Xdot[4].plot(alphas .* 180/pi, Xs[2,:], label = nothing)
    axs_X_Xdot[4].plot(alphas .* 180/pi, Xdots[2,:], label = nothing)
    axs_X_Xdot[4].plot(alphas .* 180/pi, Xdots_check[2,:], label = nothing, "--")
    axs_X_Xdot[5].plot(alphas .* 180/pi, Xs[4,:], label = nothing)
    axs_X_Xdot[5].plot(alphas .* 180/pi, Xdots[4,:], label = nothing)
    axs_X_Xdot[5].plot(alphas .* 180/pi, Xdots_check[4,:], label = nothing, "--")
    axs_X_Xdot[6].plot(alphas .* 180/pi, Xs[6,:] .* 180/pi, label = nothing)
    axs_X_Xdot[6].plot(alphas .* 180/pi, Xdots[6,:] .* 180/pi, label = nothing)
    axs_X_Xdot[6].plot(alphas .* 180/pi, Xdots_check[6,:] .* 180/pi, label = nothing, "--")

    fig_X_Xdot.legend(loc="upper left", bbox_to_anchor=(0.73,0.97))
    fig_X_Xdot.suptitle(MT.LS.L"\gamma = " * "$(gamma * 180/pi)" * MT.LS.L"^\circ")
    fig_X_Xdot.tight_layout()
    fig_X_Xdot.set_size_inches(12, 7, forward=true)

    if save_figures
        fig_X_Xdot.savefig(joinpath(MT.topdirectory, "data", "plots", MT.TODAY, fig_name))
    end

    return Xs, Xdots, Xdots_check
end

Xs, Xdots, Xdots_check = test_xdot_cl_alpha(crc3_model; gamma, fig_name = "test_dynamics_$(nwings)_wing_gamma_$(Int(round(gamma; digits = 0))).png")
_ = test_xdot_cl_alpha(crc3_model;
    alphas = range(-40.0, stop=55.0, length=26) .* pi/180,
    gamma = 0 * pi/180,
    vinf = 1.0,
    rotors_on = true,
    U = ones(length(crc3_model.aircraft.rotor_system.index)) * 3000 * 2*pi/60.0,
    save_figures = true, fig_name = "test_dynamics_$(nwings)_wing_gamma_$(Int(round(gamma; digits = 0)))_rotors_on.png", clear_figure = false, label_tag = " w/ rotors"
)
