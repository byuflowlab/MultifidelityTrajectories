#####
##### Model
#####
crc3, parameters, environment, alphas, stepsymbol, cl_stall, nx, nu = MT.build_crc3(; nwings)
crc3_model = MT.AircraftDynamics(crc3, parameters, environment, stepsymbol, cl_stall, nx, nu)

function u_2_thrust(model, U, alpha, gamma, vinf)
    # get state
    X = MT.steady_climb_state(alpha, gamma, vinf)

    # get forces/moments; aerodynamic force should be limited to thrust and drag
    force, M_y, power, rotor_force, rotor_moment = MT.verbose_aerodynamics!(model, X, U; rotors_on = true)
    println("Sherlock!\n\trotor_force = $rotor_force\n\tforce = $force")
    return force, M_y, power, rotor_force, rotor_moment
end

function test_thrust(model;
    npoints = 30,
    Us = range(0, stop=5000 .* pi/30, length=npoints) * [1,1]',
    vinf = 2.0,
    alpha = 0.0,
    gamma = pi/2
)
    forces = zeros(2,npoints)
    M_ys = zeros(npoints)
    rotor_forces = zeros(2,npoints)
    rotor_moments = zeros(npoints)

    for i in 1:size(Us)[1]
        force, M_y, power, rotor_force, rotor_moment = u_2_thrust(model, Us[i,:], alpha, gamma, vinf)
        println("\trotor_force = $rotor_force\n")
        forces[:,i] .= force
        M_ys[i] = M_y
        rotor_forces[:,i] .= rotor_force ./ 2 # just 1 per rotor
        rotor_moments[i] = rotor_moment
    end
    println("\n\nSherlock!\nrotor_forces = $rotor_forces")

    rho = model.environment.ρ
    ns = Us[:,1] ./ (2*pi)
    d = crc3_model.aircraft.rotor_system.rotors[1].Rtip * 2
    CTs_x = rotor_forces[1,:] ./ (rho * ns .^2 * d^4)
    CTs_y = rotor_forces[2,:] ./ (rho * ns .^2 * d^4)
    Js = vinf ./ (ns .* d)
    println("\nSherlock!\n\tCTs_x = $CTs_x")

    return Js, CTs_x, CTs_y, rotor_moments
end

Js, CTs_x, CTs_y, rotor_moments = test_thrust(crc3_model)

clear_figure = true
fig_thrust = PP.figure("test_thrust")
if clear_figure; fig_thrust.clear(); end
fig_thrust.add_subplot(2,1,1, ylabel = MT.LS.L"C_{T,x} [N]")
fig_thrust.add_subplot(2,1,2, ylabel = MT.LS.L"C_{T,y}", xlabel = MT.LS.L"J")
axs_thrust = fig_thrust.get_axes()

fig_moment = PP.figure("test_rotor_moment")
if clear_figure; fig_moment.clear(); end
fig_moment.add_subplot(1,1,1, ylabel = MT.LS.L"M_y [N-m]", xlabel = MT.LS.L"J")
axs_moment = fig_moment.get_axes()
# axs_thrust[2].plot(alphas .* 180/pi, CFs[1,:])
axs_thrust[1].plot(Js, CTs_x)#, label=MT.LS.L"CT_x")
axs_thrust[2].plot(Js, CTs_y)#, label=MT.LS.L"CT_y")
axs_moment[1].plot(Js, rotor_moments)

fig_thrust.savefig(joinpath(MT.topdirectory, "data", "plots", MT.TODAY, "test_thrust_forces.png"))
fig_moment.savefig(joinpath(MT.topdirectory, "data", "plots", MT.TODAY, "test_thrust_moments.png"))
