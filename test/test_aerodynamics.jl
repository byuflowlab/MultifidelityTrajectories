#####
##### Model
#####
crc3, parameters, environment, alphas, stepsymbol, cl_stall, nx, nu = MT.build_crc3(; nwings)
crc3_model = MT.AircraftDynamics(crc3, parameters, environment, stepsymbol, cl_stall, nx, nu)

#####
##### re-nondimensionalize
#####
function alpha_2_CF_CM(crc3_model, alpha; vinf = 1.0, rotors_on = false, U = [0.0, 0.0]* 2*pi/60.0, frame = MT.Wind())
    # build state
    X = MT.level_flight_state(alpha, vinf)

    # get forces/moments
    force, M_y = MT.aerodynamics!(crc3_model, X, U; rotors_on)

    # get normalization parameters
    Sref = crc3_model.aircraft.wing_system.system.reference[1].S
    b = crc3_model.aircraft.wing_system.system.reference[1].b
    c = crc3_model.aircraft.wing_system.system.reference[1].c
    rho = crc3_model.environment.ρ
    qinf = 0.5 * rho * vinf^2

    # get force/moment coefficients
    CF = force ./ qinf ./ Sref
    CM = M_y / qinf / Sref / c

    return CF, CM
end

#####
##### test!
#####
function test_cl_alpha(crc3_model; alphas = range(-20.0, stop=25.0, length=26) .* pi/180, vinf = 1.0, rotors_on = false, U = [0.0, 0.0], frame = MT.Wind(), save_figures = true, clear_figure = true)
    # initialize vars
    CFs = zeros(2,length(alphas))
    CMs = zeros(length(alphas))

    # build lift curves
    for (i,alpha) in enumerate(alphas)
        CF, CM = alpha_2_CF_CM(crc3_model, alpha; vinf, rotors_on, U, frame)
        CFs[:,i] .= CF
        CMs[i] = CM
    end

    fig_cf = PP.figure("test_CF")
    if clear_figure; fig_cf.clear(); end
    fig_cf.add_subplot(2, 1, 1, ylabel = MT.LS.L"C_L")
    fig_cf.add_subplot(2, 1, 2, ylabel = MT.LS.L"C_D", xlabel = MT.LS.L"\alpha[^\circ]")
    axs_cf = fig_cf.get_axes()
    axs_cf[2].plot(alphas .* 180/pi, CFs[1,:])
    axs_cf[1].plot(alphas .* 180/pi, CFs[2,:])

    fig_cm = PP.figure("test_CM")
    if clear_figure; fig_cm.clear(); end
    fig_cm.add_subplot(1,1,1, ylabel = MT.LS.L"C_M", xlabel = MT.LS.L"\alpha[^\circ]")
    axs_cm = fig_cm.get_axes()[1]
    axs_cm.plot(alphas .* 180/pi, CMs)

    if save_figures
        fig_cf.savefig(joinpath(MT.topdirectory, "data", "plots", MT.TODAY, "test_aerodynamics_cf_$(nwings)_wing.png"))
        fig_cm.savefig(joinpath(MT.topdirectory, "data", "plots", MT.TODAY, "test_aerodynamics_cm_$(nwings)_wing.png"))
    end

    return alphas, CFs, CMs
end

_ = test_cl_alpha(crc3_model; U=[0.0])
