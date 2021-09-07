module MultifidelityTrajectories
# this directory
const topdirectory = normpath(joinpath(@__DIR__, ".."))

# import AircraftSystems
include(joinpath(topdirectory, "..", "AircraftSystems", "src", "AircraftSystems.jl"))
const AS = AircraftSystems

# import Altro
const Dates = AS.Dates
# import Dates
# const RD = Altro.RobotDynamics
const ST = AS.StaticArrays
# const TO = Altro.TrajectoryOptimization
# using LinearAlgebra, StaticArrays
const LA = AS.LinearAlgebra
const LS = AS.LaTeXStrings
const PP = AS.PyPlot
import SNOW

# set constants
const TODAY = replace(string(Dates.today()),"-" => "")

# include files
include("CRC3_geometry.jl")
include("CRC3_system.jl")
include("coordinate_frames.jl")
include("CRC3_dynamics.jl")
include("CRC3_aerodynamics.jl")
include("utilities.jl")
include("integrate_dynamics.jl")
include("trim.jl")
include("optimization.jl")
include("plots.jl")

end # module
