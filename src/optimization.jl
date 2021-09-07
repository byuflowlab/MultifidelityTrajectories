abstract type Shooting end
struct shooting end
abstract type Collocation end
struct collocation end

function optimization_template(model, u0, xf, x0, ::Shooting;
    quadrature = rk4(),
    nt=10,
    ts=range(0.0, stop=1.0, length=nt),
    lu = zeros(length(ts) * model.nu),  # lower bounds on u
    uu = ones(length(ts) * model.nu) .* 10000 * pi/30,  # upper bounds on u
    objective = (x) -> x[end],
    solver = SNOW.IPOPT()
)
    xs = zeros(length(ts) * model.nx)
    function objective!(g, u)
        integrate_dyanmics!(xs, model, x0, u, ts, quadrature)
        # x0 constraint
        g[1:model.nx] .= xs[1:model.nx] .- x0
        # xf constraint
        g[model.nx+1:2*model.nx] .= xs[end-model.nx+1:end] .- xf

        return objective(xs)
    end

    ng = 5*model.nx  # number of constraints
    lg = -Inf*ones(ng)  # lower bounds on g
    lg[1:2*model.nx] .= 0.0 # equality constraints on start and end states
    ug = zeros(ng)  # upper bounds on g
    options = SNOW.Options(solver=solver)  # choosing IPOPT solver

    return ts, xs, objective!, u0, ng, lu, uu, lg, ug, options
end

function optimization_template(model, u0, xf, ::Shooting;
    quadrature = rk4(),
    nt=10,
    ts=range(0.0, stop=1.0, length=nt),
    lu = zeros(length(ts) * model.nu),  # lower bounds on u
    uu = ones(length(ts) * model.nu) .* 10000 * pi/30,  # upper bounds on u
    objective = (x) -> x[end],
    solver = SNOW.IPOPT()
)
    x0 = trim(model)
end

# uopt, fopt, info = SNOW.minimize(objective!, u0, ng, lu, uu, lg, ug, options)
# println("ustar = ", uopt)
# println("fstar = ", fopt)
# println("info = ", info)
