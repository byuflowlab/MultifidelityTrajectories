import MultifidelityTrajectories
const MT = MultifidelityTrajectories

##### TEST #####

#####
##### Model
#####
params = MT.build_crc3()
crc3, parameters, environment, alphas, stepsymbol, cl_stall = MT.build_crc3()
crc3_model = MT.AircraftDynamics(crc3, parameters, environment, alphas, stepsymbol, cl_stall)

#####
##### State
#####
x = y = 0.0
vinf = 1.0
alpha = 5.0 * pi/180
vx = vinf * cos(alpha)
vy = -vinf * sin(alpha)
theta = alpha
theta_dot = 0.0
X = [x, y, vx, vy, theta, theta_dot]
U = [2000, 3000] * 2*pi/60.0

Xdot = MT.RD.dynamics(crc3_model, X, U)
