crc3, parameters, environment, alphas, stepsymbol, cl_stall = MT.build_crc3(; nwings)
crc3_model = MT.AircraftDynamics(crc3, parameters, environment, stepsymbol, cl_stall)
ts = range(0, stop=0.1, length=20)
us = ones(length(ts)) .* 3000.0 * 2  * pi / 60
alpha = 5.0 * pi/180
gamma = 3.0 * pi/180
vinf = 5.0
x0 = MT.steady_climb_state(alpha, gamma, vinf)
nx = length(x0)

xs = MT.integrate_dynamics(crc3_model, x0, us, ts)
MT.plot_states(ts, xs, nx)
