
"""
    build_crc3()

Build's required objects to model the CRC-3.

# Returns

- `crc3`: `AircraftSystems.Aircraft` object
- `default_parameters`: `<: AircraftSystems.Parameters` struct
- `actions`: list of functions required to solve the aerodynamics of the CRC-3
- `alphas`: a vector containing one element to hold the angle of attack
- `stepsymbol`: required input to solve the aerodynamis using `AircraftSystems`

"""
function build_crc3()

    #####
    ##### Define wings
    #####
    wing_lexs = [
        [0.0, 0.0], # lower wing
        [0.0, 0.0]  # upper wing
    ]
    wing_leys = [
        [0.0, wing_b/2],
        [0.0, wing_b/2]
    ]
    wing_lezs = [
        [0.0, 0.0],
        [wing_spacing, wing_spacing]
    ]
    wing_chords = [
        [wing_c, wing_c],
        [wing_c, wing_c]
    ]
    wing_thetas = [
        [0.0, 0.0],
        [0.0, 0.0]
    ]
    wing_phis = [
        [0.0, 0.0],
        [0.0, 0.0]
    ]
    wing_nspanwisepanels = [40, 40]
    wing_nchordwisepanels = [1, 1]

    #####
    ##### Define rotors
    #####
    alphas = [0.0]
    omegas = ones(2) .* 2000.0 * pi/30 # rad/s
    nblades = [2]
    rtip = [rotor_R]
    radii = fill(rotor_radii, length(nblades))
    rhub = fill(rotor_rhub, length(nblades))
    chords = fill(rotor_chord, length(nblades))

    # add collective
    twists = fill(rotor_twist, length(nblades)) # radians
    twist_extra = 0.0 * pi/180
    for (i, twist) in enumerate(twists)
        twist .+= twist_extra
    end

    # extra geometry
    index = [1,1]
    rotor_positions = [
        [-rotor_x, rotor_y, 0.0],
        [-rotor_x, rotor_y, wing_spacing]
    ]
    rotor_orientation = fill([-1.0 * cos(rotor_cant), -sin(rotor_cant), 0.0], length(index))
    spindirections = fill(true, length(index))

    #####
    ##### define airfoils
    #####
    # contours
    contourdirectory = joinpath(topdirectory, "data", "airfoil", "contour")
    contourfilenames = fill("e212",length(radii[1]))
    airfoilcontours = fill(joinpath.(contourdirectory, contourfilenames .* ".dat"),length(nblades))

    # polars
    polar_directory = joinpath(topdirectory, "data", "airfoil", "polar")
    airfoil_name = "crc3_rR06"
    polar_paths = [joinpath(topdirectory, "data", "airfoil", "polar", "contour66.33_polar_5000.0_0.62.csv")]
    Res = [5000]
    Ms = [0.1]
    viternaextrapolation = false
    rotationcorrection = false
    AS.rename_polars(polar_paths, Res, Ms, airfoil_name, viternaextrapolation, rotationcorrection, polar_directory; force = true)
    airfoilnames = fill(fill(airfoil_name,length(radii[1])), length(nblades))

    plotstepi = 1:length(alphas)

    Res = [5e3]
    Machs = [0.1]

    Res_list = fill(fill(Res, length(radii[1])), length(nblades))
    Ms_list = fill(fill(Machs, length(radii[1])), length(nblades))

    vinfs = [1.0] # dummy v

    #####
    ##### build system
    #####
    crc3, _, _, _, _, _, _, alphas, stepsymbol = AS.vlm_bem_template(vinfs, plotstepi, alphas,
        wing_lexs, wing_leys, wing_lezs, wing_chords, wing_thetas, wing_phis, wing_nspanwisepanels, wing_nchordwisepanels,
        omegas, nblades, rhub, rtip, radii, chords, twists,
        airfoilcontours, airfoilnames, index,
        rotor_positions, rotor_orientation, spindirections, Res_list, Ms_list;
            Vref = 1.0, symmetric = true, iref = 1, staticmargin = 0.10, liftingline_x_over_c = 0.25,
            wing_cambers = [fill((xc) -> 0, length(lex)) for lex in wing_lexs],
            spanwisespacings = fill(AS.VL.Uniform(), length(wing_lexs)),
            chordwisespacings = fill(AS.VL.Uniform(), length(wing_lexs)),
            wakedevelopementfactor = 1.0, # fully developed by default
            swirlrecoveryfactor = 0.5, # as described in Veldhuis' paper
            surfacenames = ["lower wing", "upper wing"],
            rotornames = ["lower rotor", "upper rotor"],
            mass = mass_total,
            inertia_x = 0.0,
            inertia_y = JInertia,
            inertia_z = 0.0,
            plotdirectory = joinpath(topdirectory, "data","plots",TODAY),
            plotbasename = "crc3",
            plotextension = ".png",
            stepsymbol = LS.L"\alpha ",
            polardirectory = polar_directory
    )

    #####
    ##### parameters
    #####
    solve_rotor_nondimensional_params = AS.solve_rotor_nondimensional(crc3, alphas)
    solve_rotor_wake_params = AS.solve_rotor_wake(crc3, alphas)
    solve_rotor_wake_params = solve_rotor_wake_params[[1,4,5,6,7,8]]
    solve_wing_CF_CM_params = AS.solve_wing_CF_CM(crc3, alphas)[2:end]
    lift_moment_distribution_params = AS.lift_moment_distribution(crc3, alphas)

    nparams = length(solve_rotor_nondimensional_params) + length(solve_rotor_wake_params) + length(solve_wing_CF_CM_params) + length(lift_moment_distribution_params)
    println("nparams = $nparams")
    all_params = vcat(solve_rotor_nondimensional_params, solve_rotor_wake_params, solve_wing_CF_CM_params, lift_moment_distribution_params)
    println("types = $(typeof.(all_params))")

    parameters = AS.ParamsSolveVLMBEM(solve_rotor_nondimensional_params..., solve_rotor_wake_params..., solve_wing_CF_CM_params..., lift_moment_distribution_params...)

    #####
    ##### Environment
    #####
    environment = AS.Environment()

    #####
    ##### Stall c_l
    #####
    cl_stall = 1.3

    return crc3, parameters, environment, alphas, stepsymbol, cl_stall
end
