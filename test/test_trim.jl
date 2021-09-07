import Snopt

#####
##### Model
#####

println("===== Building Model =====")
crc3, parameters, environment, alphas, stepsymbol, cl_stall, nx, nu = MT.build_crc3(; nwings)
crc3_model = MT.AircraftDynamics(crc3, parameters, environment, stepsymbol, cl_stall, nx, nu)

#####
##### trim
#####

# xu_xdot = [
#     1. px,
#     2. py,
#     3. vx,
#     4. vy,
#     5. theta,
#     6. theta_dot,
#     7. energy,
#     8. px_dot,
#     9. py_dot,
#     10. vx_dot,
#     11. vy_dot,
#     12. theta_dot,
#     13. theta_dot_dot,
#     14. power
# ]

x_trim_guess = [
    0.0, # x
    0.0, # y
    0.0, # vx
    0.0, # vy
    pi/2, # theta
    0.0, # theta_dot
    0 # energy
]

u_trim_guess = ones(nwings) * 10000 .* pi/30

xu_trim_guess = vcat(x_trim_guess, u_trim_guess)

x_xdot_trimmed = [
    # 0.0, # px
    # 0.0, # py
    # 0.0, # vx
    # 0.0, # vy
    pi/2, # theta
    # 0.0, # theta_dot
    # 0.0, # energy
    # 0.0, # px_dot
    # 0.0, # py_dot
    0.0, # vx_dot
    0.0, # vy_dot
    # 0.0, # theta_dot
    0.0, # theta_dot_dot
    # 0.0, # power
]

i_x_xdot_trimmed = [5,10,11,13]
# i_xu_xdot_trimmed = [3, 4, 6, 7, 8, 9, 10, 11, 12, 13]

lxu = [
    0.0,
    0.0,
    -10.0,
    -10.0,
    -pi,
    -10,
    0.0,
    0.0,
    0.0
]

uxu = [
    0.0,
    0.0,
    10.0,
    10.0,
    pi,
    10.0,
    10.0,
    10000.0 * pi/30,
    10000.0 * pi/30
]

println("===== Begin Optimization =====")

xopt, fopt, info = MT.trim(crc3_model, xu_trim_guess, x_xdot_trimmed, i_x_xdot_trimmed, lxu, uxu;
    derivatives = MT.SNOW.CentralFD(),
    # derivatives = MT.SNOW.ComplexStep(),
    # solver = MT.SNOW.SNOPT(),
    solver = MT.SNOW.IPOPT(),
    objective = (x, xdot, u) -> MT.LA.norm(x[3:4]) / 1000
)

println("xstar = ", xopt)
println("fstar = ", fopt)
println("info = ", info)
