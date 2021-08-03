module MultifidelityTrajectories
# this directory
const topdirectory = normpath(joinpath(@__DIR__, ".."))

# import AircraftSystems
include(joinpath(topdirectory, "..", "AircraftSystems", "src", "AircraftSystems.jl"))
const AS = AircraftSystems

import Altro
const Dates = AS.Dates
# import Dates
const RD = Altro.RobotDynamics
const SArr = Altro.StaticArrays
const TO = Altro.TrajectoryOptimization
# using LinearAlgebra, StaticArrays
const LA = AS.LinearAlgebra
const LS = AS.LaTeXStrings

# set constants
const TODAY = replace(string(Dates.today()),"-" => "")

# include files
include("CRC3_geometry.jl")
include("CRC3_system.jl")
include("CRC3_dynamics.jl")

end # module
