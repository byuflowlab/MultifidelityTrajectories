"""
    trim(model, xu_trim_guess, xu_xdot_trimmed, i_xu_xdot_trimmed, lxu, uxu;
        derivatives = SNOW.ComplexStep(),
        solver = SNOW.IPOPT(),
        objective = (g, xu) -> LA.norm(g, Inf)
    )

Obtains a trimmed state of the model.

# Arguments

* `model<:AbstractModel` - model to be trimmed
* `xu_trim_guess::Vector{Float64}` - guess as to the concatenated state and controls when trimmed
* `xdot_trimmed::Vector{Float64}` - values taken on by any number of the time derivative of the state when trimmed
* `i_xdot_trimmed::Vector{Int}` - indices of the time derivative of the state corresponding to `xdot_trimmed` (these must be the same length)
* `lxu::Vector{Float64}` - lower bounds of the concatenated state and control vector
* `uxu::Vector{Float64}` - upper bounds of the concatenated state and control vector

# Keyword Arguments

* `derivatives = SNOW.ComplexStep()`
* `solver = SNOW.IPOPT()`
* `objective = (g, xu) -> LA.norm(g, Inf)`

"""
function trim(model, xu_trim_guess, x_xdot_trimmed, i_x_xdot_trimmed, lxu, uxu;
    derivatives = SNOW.ComplexStep(),
    solver = SNOW.IPOPT(),
    objective = (x, xdot, u) -> xdot[7]
)
    @assert length(x_xdot_trimmed) == length(i_x_xdot_trimmed) "length of xu_xdot_trimmed and i_xu_xdot_trimmed inconsistent"
    @assert length(xu_trim_guess) == model.nx + model.nu "length of xu_trim_guess must be the sum of states and controls nx + nu"

    # prepare optimizer options
    options = SNOW.Options(derivatives = derivatives, solver = solver)
    ng = length(x_xdot_trimmed) * 2
    lg = fill(-Inf, ng)
    # lg = zeros(ng) # equality constraint
    ug = zeros(ng)
    xdot = zeros(model.nx)

    function objective!(g, xu)
        x = xu[1:model.nx]
        u = xu[model.nx+1:end]
        xdot .= dynamics(model, x, u)
        x_xdot = vcat(x,xdot)
        g[1:Int(ng/2)] .= x_xdot[i_x_xdot_trimmed] .- x_xdot_trimmed
        g[1+Int(ng/2):end] .= x_xdot_trimmed .- x_xdot[i_x_xdot_trimmed]

        return objective(x, xdot, u)
    end

    # println("Sherlock!\n\txu_trim_guess = $xu_trim_guess\n\tng = $ng\n\ttlxu = $lxu\n\tuxu = $uxu\n\tlg = $lg\n\tug = $ug\n\toptions = $options")
    # println("\nTest ===:")
    # for i=1:length(xu_trim_guess)
    #     println(xu_trim_guess[i] .=== xu_trim_guess)
    # end
    xu_opt, f_opt, info = SNOW.minimize(objective!, xu_trim_guess, ng, lxu, uxu, lg, ug, options)

    return xu_opt, f_opt, info
end
