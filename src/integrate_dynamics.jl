abstract type Quadrature end
abstract type Euler <: Quadrature end
struct euler <: Euler end
abstract type RK4 <: Quadrature end
struct rk4 <: RK4 end

function integrate_dynamics(model, x0, us, ts, quadrature = euler())
    xs = Array{typeof(x0[1]),1}(undef, length(x0) * length(ts))
    integrate_dynamics!(xs, model, x0, us, ts, quadrature)
    return xs
end

"""
    integrate_dyanmics!(xs, model, x0, us, ts, quadrature)

Integrates the trajectory provided control inputs. Modifies `xs` in place.

# Inputs

* `xs::Vector{Float64}` - preallocated memory consisting of concatenated state vectors to be populated with simulation states for each element of `ts`
* `model<:AbstractModel` - model to be simulated
* `x0::Vector{Float64}` - initial state
* `us::Vector{Float64}` - concatenated control inputs prescribed for the simulation for each element of `ts`
* `ts::Vector{Float64}` - vector of times at which dynamics are evaluated
* `quadrature<:Quadrature` - what kind of quadrature is to be used for integration; `euler()` and `rk4()` currently supported

"""
function integrate_dynamics!(xs, model, x0, us, ts, quadrature::Euler)
    nx = model.nx # number of states
    @assert nx == length(x0) "length of x0 and model.nx inconsistent"
    nu = model.nu # number of controls
    @assert nu == Int(length(us) / length(ts)) "length of us and length of ts and model.nu inconsistent"
    @assert length(xs) == length(x0) * length(ts) "length of states x not consistent with timesteps and x0"

    # initial knot is the initial guess
    xs[1:nx] .= x0

    # loop through time
    for i=1:length(ts)-1
        Δt = ts[i+1] - ts[i]
        x = xs[(i-1) * nx + 1 : i * nx]
        u = us[(i-1) * nu + 1 : i * nu]
        xdot = dynamics(model, x, u)
        xs[i * nx + 1 : (i+1) * nx] .= x + xdot * Δt
    end

    return nothing
end

function integrate_dyanmics!(xs, model, x0, us, ts, quadrature::RK4)
    nx = model.nx # number of states
    @assert nx == length(x0) "length of x0 and model.nx inconsistent"
    nu = model.nu # number of controls
    @assert nu == Int(length(us) / length(ts)) "length of us and length of ts and model.nu inconsistent"
    @assert length(xs) == length(x0) * length(ts) "length of states x not consistent with timesteps and x0"

    # initial knot is the initial guess
    xs[1:nx] .= x0

    # loop through time
    for i=1:length(ts)-1
        Δt = ts[i+1] - ts[i]
        x = xs[(i-1) * nx + 1 : i * nx]
        u = us[(i-1) * nu + 1 : i * nu]
        u_next = us[i * nu + 1 : (i+1) * nu]
        F1 = dynamics(model, x, u)
        F2 = dynamics(model, x + .5 * Δt .* F1, u)
        F3 = dynamics(model, x + .5 * Δt .* F2, u_next)
        F4 = dynamics(model, x +  Δt .* F3, u_next)
        xs[i * nx + 1 : (i+1) * nx] .= x + Δt/6 .* (F1 + 2 .* F2 + 2 .* F3 + F4)
    end

    return nothing
end
