#####
##### dynamics model
#####
struct AircraftDynamics{TF, TAF}
    aircraft::AS.Aircraft{TF,TF,TF,TF,TAF}
    parameters::AS.ParamsSolveVLMBEM{TF}
    environment::AS.Environment{TF}
    stepsymbol::Union{String, LS.LaTeXString}
    cl_stall::TF
    nx::Int64
    nu::Int64
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
    * `X[7]` energy expended [J]

* `U::Vector{Float64}` vector of control inputs in the following order

    * `U[1]` radians/second of the lower rotor
    * `U[2]` radians/second of the upper rotor

# Returns

* `X_dot::Vector{Float64}` vector of time derivatives of states as follows:

    * `X[1]` x_dot; x component of velocity in the global frame
    * `X[2]` y_dot; y component of velocity in the global frame
    * `X[3]` vx_dot; x component of acceleration in the body frame
    * `X[4]` vy_dot; y component of acceleration in the body frame
    * `X[5]` theta_dot; time rate of change of theta
    * `X[6]` theta_dot_dot; time acceleration of theta
    * `X[7]` instantaneous power [W]

"""
function dynamics(model, X, U; rotors_on = true)
    # unpack parameters
    aircraft = model.aircraft
    environment = model.environment # assumed constant
    g = environment.g
    mass = aircraft.inertia_system.mass
    inertia = aircraft.inertia_system.inertia_y

    # calculate resultants in the body frame
    force, M_y, power = aerodynamics!(model, X, U; rotors_on, frame = Body())

    # # calculate power
    # Qs = parameters.Qs[:,1]
    # Ps = Qs .* U
    # P = 2 * sum(Ps)

    return dynamics(X, mass, inertia, g, force, M_y, power)
end

@inline function dynamics(X, mass, inertia, g, force, moment, power)
    theta = X[5]

    s, c = sincos(theta)
    R = ST.@SArray [c -s; s c] # rotation from body to global frame
    v = ST.@SVector [X[3], X[4]]
    theta_dot = X[6]

    p_dot = R * v
    v_dot = ( -g * R' * (ST.@SVector [0,1])  # apply gravity in inertial -z direction
            + force / mass                  # aerodynamic forces
            - theta_dot .* (ST.@SArray [0. -1.; 1. 0.]) * v) # coriolis forces
    theta_dot_dot = moment / inertia

    return ST.@SVector [p_dot[1], p_dot[2], v_dot[1], v_dot[2], theta_dot, theta_dot_dot, power]
end
