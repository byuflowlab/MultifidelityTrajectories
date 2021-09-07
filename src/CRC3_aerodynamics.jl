"""
    body_thrust(Ts, orientations; rotors_on = true)

Returns thrust in the body frame in 3-D.
"""
@inline function body_thrust(Ts, orientations; rotors_on = true)
    Ts_vec = Vector{ST.SArray{Tuple{3},Float64,1,3}}(undef,length(Ts)) # 2-D
    # Ts_vec = Vector{ST.SArray{Tuple{3},Float64,1,3}}(undef,length(Ts)) # 3-D
    if rotors_on
        for (i,T) in enumerate(Ts)
            # T_vec = ST.@SVector [orientations[i][1] * T, orientations[i][3]] # 2-D
            Ts_vec[i] = ST.@SVector [orientations[i][1] * T, orientations[i][2] * T, orientations[i][3] * T] # 3-D
        end
    else
        for (i,T) in enumerate(Ts)
            Ts_vec[i] = ST.@SVector [0.0, 0.0, 0.0] # 2-D
            # Ts_vec[i] = ST.@SVector [orientations[i][1] * T, orientations[i][2] * T, orientations[i][3]] 3-D
        end
    end

    return Ts_vec
end

function state_2_freestream(model::AircraftDynamics, X)
    # calculate freestream
    vx, vy, theta, theta_dot = X[3:6]
    airspeed = [vx, 0.0, -vy] .+ [1e-12, 0.0, 0.0] # in body frame
    vinf = LA.norm(airspeed)
    alpha = asin(airspeed[3] / vinf)
    Omega = ST.@SVector [0.0, theta_dot, 0.0]
    beta = asin(airspeed[2] / vinf)
    freestream = AS.Freestream(vinf, alpha, beta, Omega)

    return freestream
end

function rotor_force(aircraft, parameters, X, freestream, frame; rotors_on = true)
    # extract rotor performance
    Ts = parameters.Ts[:,1]
    orientations = aircraft.rotor_system.orientations
    Ts_vec = body_thrust(Ts, orientations; rotors_on) # in the desired frame, in 3-D

    # calculate forces and moments caused by the rotor in the specified frame
    theta = X[5]
    alpha = freestream.alpha
    Frotor = R_from_body(frame, alpha, theta) * sum(Ts_vec)[[1,3]] * 2 # *2 to mirror
    positions = aircraft.rotor_system.positions
    @assert length(positions) == length(Ts_vec) "Length of rotor positions ($(length(positions))) and thrust vectors ($(length(Ts_vec))) inconsistent."
    # y-axis moment doesn't change between frames
    Mrotor = sum([LA.cross(positions[i] - aircraft.wing_system.system.reference[1].r, Ts_vec[i])[2] for i in 1:length(Ts_vec) ]) # (Ts_vec[2][1] - Ts_vec[1][1]) * (aircraft.wing_system.system.surfaces[2][1].rcp[3] - aircraft.wing_system.system.surfaces[1][1].rcp[3])

    return Frotor, Mrotor
end

"""
    aerodynamics!(model::AircraftDynamics, X, U; rotors_on = true, frame = body::Frame)

Returns x and z forces and y moment in the specified frame.

# Arguments

* `model::AircraftDynamics`
* `X::Vector{Float64}`
* `U::Vector{Float64}`

# Keyword Arguments

* `rotors_on::Bool = true`
* `frame::Frame = Body()`

# Returns

* force
* moment
* power

"""
function aerodynamics!(model::AircraftDynamics, X, U; rotors_on = true, frame::Frame = Body())
    # unpack inputs to AircraftSystems
    aircraft = model.aircraft
    parameters = model.parameters
    environment = model.environment # assumed constant
    stepi = 1 # since there is only 1 step
    stepsymbol = model.stepsymbol

    # unpack states and controls
    x, y, vx, vy, theta, theta_dot, energy = X
    # println("Sherlock!\n\tomegas = $(parameters.omegas)\n\tU = $U")
    parameters.omegas[:,1] .= U # rad/s of the rotors

    # unpack other parameters
    cl_stall = model.cl_stall
    Sref = aircraft.wing_system.system.reference[1].S
    b = aircraft.wing_system.system.reference[1].b
    c = aircraft.wing_system.system.reference[1].c
    g = environment.g

    freestream = state_2_freestream(model, X)

    # solve aerodynamics
    AS.solve_vlm_bem(aircraft, parameters, freestream, environment, (freestream.alpha), stepi, stepsymbol) # solves in the wind axis frame

    # extract aerodynamic forces and moments in wind frame
    CF = ST.@SVector [parameters.CFs[1,1], parameters.CFs[2,1], parameters.CFs[3,1]]
    CM = ST.@SVector [parameters.CMs[1,1], parameters.CMs[2,1], parameters.CMs[3,1]]

    # redimensionalize
    ρ = environment.ρ
    qinf = 0.5 * ρ * freestream.vinf^2

    F = CF .* qinf .* Sref
    M = CM .* (ST.@SVector [b, c, b]) * qinf * Sref

    Qs = parameters.Qs[:,1]
    Ps = Qs .* U
    P = 2 * sum(Ps)

    Frotor, Mrotor = rotor_force(aircraft, parameters, X, freestream, frame; rotors_on)

    # calculate resultants in the specified frame
    force = R_from_wind(frame, freestream.alpha, theta) * F[[1,3]] + Frotor # add thrust
    # force = (ST.@SVector F[[1,3]]) + Frotor # get force vector in x,z plane
    M_y = M[2] + Mrotor

    return force, M_y, P
end

"""
    verbose_aerodynamics!(model::AircraftDynamics, X, U; rotors_on = true, frame::Frame = Body())

Aerodynamics function returns rotor thrust and moment as well.

# Arguments

* `model::AircraftDynamics`
* `X::Vector{Float64}`
* `U::Vector{Float64}`

# Keyword Arguments

* `rotors_on::Bool = true`
* `frame::Frame = Body()`

# Returns

* `force::Vector{Float64}`
* `moment::Float64`
* `power::Float64`
* `rotor_force::Vector{Float64}`
* `rotor_moment::Float64`

"""
# function verbose_aerodynamics!(model::AircraftDynamics, X, U; rotors_on = true, frame::Frame = Body())
#     freestream = state_2_freestream(model, X)
#     total_force, total_moment, power = aerodynamics!(model, X, U; rotors_on, frame)
#     rotor_f, rotor_m = rotor_force(model.aircraft, model.parameters, X, freestream, frame; rotors_on)

#     return total_force, total_moment, power, rotor_f, rotor_m
# end
function verbose_aerodynamics!(model::AircraftDynamics, X, U; rotors_on = true, frame::Frame = Body())
    # unpack inputs to AircraftSystems
    aircraft = model.aircraft
    parameters = model.parameters
    environment = model.environment # assumed constant
    stepi = 1 # since there is only 1 step
    stepsymbol = model.stepsymbol

    # unpack states and controls
    x, y, vx, vy, theta, theta_dot, energy = X
    # println("Sherlock!\n\tomegas = $(parameters.omegas)\n\tU = $U")
    parameters.omegas[:,1] .= U # rad/s of the rotors

    # unpack other parameters
    cl_stall = model.cl_stall
    Sref = aircraft.wing_system.system.reference[1].S
    b = aircraft.wing_system.system.reference[1].b
    c = aircraft.wing_system.system.reference[1].c
    g = environment.g

    freestream = state_2_freestream(model, X)

    # solve aerodynamics
    AS.solve_vlm_bem(aircraft, parameters, freestream, environment, (freestream.alpha), stepi, stepsymbol) # solves in the wind axis frame

    # extract aerodynamic forces and moments in wind frame
    CF = ST.@SVector [parameters.CFs[1,1], parameters.CFs[2,1], parameters.CFs[3,1]]
    CM = ST.@SVector [parameters.CMs[1,1], parameters.CMs[2,1], parameters.CMs[3,1]]

    # redimensionalize
    ρ = environment.ρ
    qinf = 0.5 * ρ * freestream.vinf^2

    F = CF .* qinf .* Sref
    M = CM .* (ST.@SVector [b, c, b]) * qinf * Sref

    Qs = parameters.Qs[:,1]
    Ps = Qs .* U
    P = 2 * sum(Ps)

    Frotor, Mrotor = rotor_force(aircraft, parameters, X, freestream, frame; rotors_on)

    # calculate resultants in the specified frame
    force = R_from_wind(frame, freestream.alpha, theta) * F[[1,3]] + Frotor # add thrust
    # force = (ST.@SVector F[[1,3]]) + Frotor # get force vector in x,z plane
    M_y = M[2] + Mrotor

    println("Sherlock!\n\tparameters.CTs = $(parameters.CTs)\n\tparameters.Ts = $(parameters.Ts)\n\tparameters.omegas = $(parameters.omegas)\n\tvinf = $(freestream.vinf)")
    println("\tC_T = $(parameters.Ts[1] ./ (ρ * (parameters.omegas[1] ./ (2*pi))^2 * (2 * aircraft.rotor_system.rotors[1].Rtip)^4))")

    return force, M_y, P, Frotor, Mrotor
end
