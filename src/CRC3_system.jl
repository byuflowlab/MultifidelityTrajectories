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
function build_crc3(; nwings = 2)

    #####
    ##### Define wings
    #####

    wing_lexs = fill([0.0, 0.0], nwings)
    wing_leys = fill([0.0, crc3_wing_b/2], nwings)
    wing_lezs = [[crc3_wing_spacing, crc3_wing_spacing] .* i for i in 0:nwings-1]
    wing_chords = fill([crc3_wing_c, crc3_wing_c], nwings)
    wing_thetas = fill([0.0, 0.0], nwings)
    wing_phis = fill([0.0, 0.0], nwings)
    wing_nspanwisepanels = fill(40, nwings)
    wing_nchordwisepanels = fill(1, nwings)

    #####
    ##### Define rotors
    #####

    alphas = [0.0]
    omegas = ones(nwings) .* 2000.0 * pi/30 # rad/s
    nblades = [2]
    rtip = [crc3_rotor_R]
    radii = fill(crc3_rotor_radii, length(nblades))
    rhub = fill(crc3_rotor_rhub, length(nblades))
    chords = fill(crc3_rotor_chord, length(nblades))

    # add collective
    twists = fill(crc3_rotor_twist, length(nblades)) # radians
    twist_extra = 0.0 * pi/180
    for (i, twist) in enumerate(twists)
        twist .+= twist_extra
    end

    # extra geometry
    index = fill(1, nwings)
    rotor_positions = [[-crc3_rotor_x, crc3_rotor_y, crc3_wing_spacing * i] for i in 0:nwings-1]
    rotor_orientations = fill([-1.0 * cos(crc3_rotor_cant), -sin(crc3_rotor_cant), 0.0], nwings)
    spin_directions = fill(true, nwings)

    #####
    ##### define airfoils
    #####
    Res = [5e3]
    Machs = [0.1]

    Res_lists = fill(Res,length(nblades))
    Ms_lists = fill(Machs,length(nblades))

    # contours
    contourdirectory = joinpath(topdirectory, "data", "airfoil", "contour")
    contour_paths = fill(joinpath(contourdirectory, "e212.dat"), length(nblades))
    # contour_paths = fill(fill(joinpath(contourdirectory, "e212.dat"), length(Res), length(Res)),length(nblades))

    # polars
    polar_directory = joinpath(topdirectory, "data", "airfoil", "polar")
    uncorrected_polar_paths = fill([joinpath(polar_directory, "crc3_rR06_Re5000_M0.1.txt")], length(nblades))
    corrected_polar_paths = fill([joinpath(polar_directory, "crc3_rR06_Re5000_M0.1_ext_rot.txt")], length(nblades))
    # uncorrected_polar_paths = fill([joinpath(polar_directory, "contour66.33_polar_5000.0_0.62.csv")], length(nblades))
    # corrected_polar_paths = fill([joinpath(polar_directory, "contour66.33_polar_5000.0_0.62_ext_rot.csv")], length(nblades))
    Res = [5000]
    Ms = [0.1]
    viterna_extrapolation = false
    rotation_correction = false

    plotstepi = 1:length(alphas)

    # Res_lists = fill(fill(Res, length(radii[1])), length(nblades))
    # Ms_lists = fill(fill(Machs, length(radii[1])), length(nblades))

    vinfs = [1.0] # dummy v

    #####
    ##### build system
    #####

    crc3, _, _, _, _, _, _, alphas, stepsymbol = AS.vlm_bem_template(vinfs, plotstepi, alphas,
        wing_lexs, wing_leys, wing_lezs, wing_chords, wing_thetas, wing_phis, wing_nspanwisepanels, wing_nchordwisepanels,
        omegas, nblades, rhub, rtip, radii, chords, twists,
        contour_paths, uncorrected_polar_paths, corrected_polar_paths, index,
        rotor_positions, rotor_orientations, spin_directions, Res_lists, Ms_lists;
            Vref = 1.0, symmetric = true, iref = 1, static_margin = 0.10, liftingline_x_over_c = 0.25,
            cambers = [fill((xc) -> 0, length(lex)) for lex in wing_lexs],
            spanwise_spacings = fill(AS.VL.Uniform(), length(wing_lexs)),
            chordwise_spacings = fill(AS.VL.Uniform(), length(wing_lexs)),
            wake_developement_factor = 1.0, # fully developed by default
            swirl_recovery_factor = 0.5, # as described in Veldhuis' paper
            surfacenames = ["default wing"],
            rotor_names = ["rotor 1"],
            mass = crc3_mass,
            inertia_x = 0.0,
            inertia_y = crc3_JInertia,
            inertia_z = 0.0,
            plot_directory = joinpath(topdirectory, "data","plots",TODAY),
            plot_base_name = "default",
            plot_extension = ".pdf",
            step_symbol = LS.L"\alpha "
    )

    #####
    ##### parameters
    #####

    solve_rotor_nondimensional_params = AS.solve_rotor_nondimensional(crc3, alphas)
    solve_rotor_wake_params = AS.solve_rotor_wake(crc3, alphas)
    solve_rotor_wake_params = solve_rotor_wake_params[[1,4,5,6,7,8]]
    solve_wing_CF_CM_params = AS.solve_wing_CF_CM(crc3, alphas)[2:end]
    lift_moment_distribution_params = AS.lift_moment_distribution(crc3, alphas)

    parameters = AS.ParamsSolveVLMBEM(solve_rotor_nondimensional_params..., solve_rotor_wake_params..., solve_wing_CF_CM_params..., lift_moment_distribution_params...)

    #####
    ##### Environment
    #####

    environment = AS.Environment()

    #####
    ##### Stall c_l
    #####

    cl_stall = 1.3

    #####
    ##### state/control
    #####

    nx = 7
    nu = 2

    return crc3, parameters, environment, alphas, stepsymbol, cl_stall, nx, nu
end
