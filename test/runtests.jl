import MultifidelityTrajectories
MT = MultifidelityTrajectories
# const MT = MultifidelityTrajectories
import PyPlot
PP = PyPlot

##### TEST #####
nwings = 2
gamma = 0.0 * pi/180

# include("test_aerodynamics.jl")
# include("test_dynamics.jl")
# include("test_integration.jl")
include("test_trim.jl")
