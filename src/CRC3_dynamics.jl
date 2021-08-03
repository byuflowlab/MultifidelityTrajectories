#####
##### dynamics model
#####
struct AircraftDynamics{TF,TS} <: RD.AbstractModel
    aircraft::AS.Aircraft
    parameters
    environment::AS.Environment
    alphas::Vector{TF}
    stepsymbol::TS
    cl_stall::TF
end

"""
    dynamics(model::Aircraft, X, U)

Describes dynamic response of an aircraft built in AircraftSystems.

# Inputs

* `model::Aircraft` inherits from `RobotDynamics.AbstractModel`
* `X::Vector{Float64}` vector of states in the following order:

    * `X[1]` x; position in the global frame
    * `X[2]` y; position in the global frame
    * `X[3]` vx; x component of velocity in the body frame
    * `X[4]` vy; y component of velocity in the body frame
    * `X[5]` theta; angle from the horizon to the aircraft's attitude
    * `X[6]` theta_dot; time rate of change of theta
    # * `X[7]` p_lower; power consumption of the lower rotor
    # * `X[8]` p_upper; power consumption of the upper rotor
    # * `X[9]` cl_max_lower; maximum local lift coefficient of the lower wing
    # * `X[10]` cl_max_upper; maximum local lift coefficient of the upper wing

* `U:Vector{Float64}` vector of control inputs in the following order

    * `U[1]` RPM of the lower rotor
    * `U[2]` RPM of the upper rotor

"""
function RD.dynamics(model::AircraftDynamics, X, U)
    # unpack states and controls
    x, y, vx, vy, theta, theta_dot = X
    u_upper, u_lower = U

    # unpack inputs to AircraftSystems
    aircraft = model.aircraft
    parameters = model.parameters
    environment = model.environment # assumed constant
    alphas = model.alphas # just a place holder for the angle of attack
    @assert length(alphas) == 1 "alphas must be a vector of length 1"
    stepi = 1 # since there is only 1 step
    stepsymbol = model.stepsymbol

    # unpack other parameters
    cl_stall = model.cl_stall
    Sref = aircraft.wingsystem.system.reference[1].S
    g = environment.g
    mass = aircraft.inertiasystem.mass
    inertia = aircraft.inertiasystem.inertia_y

    # calculate freestream
    airspeed = [vx, 0.0, -vy] .+ [1e-12, 0.0, 0.0]
    vinf = LA.norm(airspeed)
    ρ = environment.ρ
    qinf = 0.5 * ρ * vinf^2
    alphas[1] = asin(airspeed[3] / vinf)
    Omega = [0.0; theta_dot; 0.0]
    beta = 0.0
    freestream = AS.Freestream(vinf, alphas[1], beta, Omega)

    # solve aerodynamics
    println("Sherlock! aerodynamics\n\tairspeed = $airspeed\n\tU = $U\n\tJ = $(vinf / rotor_R / (U[1] /60.0))")
    AS.solve_vlm_bem(aircraft, parameters, freestream, environment, alphas, stepi, stepsymbol)

    # extract aerodynamic forces and moments
    CF = SArr.@SVector [parameters.CDs[1], parameters.CYs[1], parameters.CLs[1]]
    CM = SArr.@SVector [parameters.CMxs[1], parameters.CMys[1], parameters.CMzs[1]]

    F = CF .* qinf .* Sref
    M = CM .* (SArr.@SVector [b, c, b]) * qinf * Sref
    println("Sherlock! aerodynamic forces\n\tF = $F\n\tM = $M")

    # extract rotor performance
    Ts = parameters.Ts[:,1]
    Ts_vec = [SArr.@SVector ((T * aircraft.rotorsystem.orientations[i])[[1,3]]) for (i,T) in enumerate(Ts)]
    Qs = parameters.Qs[:,1]
    Ps = Qs .* U
    P = 2 * sum(Ps)
    println("Sherlock! rotor performance\n\tTs = $Ts\n\tTs_vec = $Ts_vec\n\tQs = $Qs\n\tPs = $Ps\n\tP = $P")

    # calculate forces and moments caused by the rotor
    Frotor = sum(Ts) * 2 # mirror
    Mrotor = (Ts[1] - Ts[2]) * (aircraft.wingsystem.system.surfaces[2][1].rcp[3] - aircraft.wingsystem.system.surfaces[1][1].rcp[3])
    println("Sherlock! rotor forces and moments\n\tFrotor = $Frotor\n\tMrotor = $Mrotor")

    # calculate resultants
    force = (SArr.@SVector F[[1,3]]) + Frotor # get force vector in x,z plane
    M_y = M[2] + Mrotor
    println("Sherlock! resultants\n\tforce = $force\n\tM_y = $M_y")

    s = sin(theta)
    c = cos(theta)
    R = SArr.@SArray [c -s; s c] # rotation from body to inertial
    v = SArr.@SVector [vx, vy]

    p_dot = R * v
    v_dot = ( -g * R' * (SArr.@SVector [0,1])  # apply gravity in inertial -z direction
            + force / mass                  # aerodynamic forces
            - theta_dot .* (SArr.@SArray [0. -1.; 1. 0.]) * v) # coriolis forces
    theta_dot_dot = M_y / inertia

    # derivatives not needed; set to zero
    p_lower_dot, p_upper_dot, cl_max_lower_dot, cl_max_upper_dot = 0.0, 0.0, 0.0, 0.0

    return [p_dot..., v_dot..., theta_dot, theta_dot_dot]
end
