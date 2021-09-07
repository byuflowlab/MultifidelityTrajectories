function index_states(xs, nx, i)
    nt = Int(length(xs) / nx)
    this_state = [xs[(j-1) * nx + i] for j in 1:nt]
    return this_state
end

function level_flight_state(alpha, vinf)
    x = y = 0.0
    vx = vinf * cos(alpha)
    vy = -vinf * sin(alpha)
    theta = alpha
    theta_dot = 0.0
    power = 0.0
    return [x, y, vx, vy, theta, theta_dot, power]
end


function steady_climb_state(alpha, gamma, vinf)
    theta = gamma + alpha
    x = y = 0.0
    v_global = ST.@SVector [vinf * cos(gamma), vinf * sin(gamma)]
    v_body = R_from_global(Body(), alpha, theta) * v_global
    theta_dot = 0.0
    energy = 0.0
    return vcat(x, y, v_body, theta, theta_dot, energy)
end
