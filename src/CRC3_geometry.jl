const crc3_rotor_R = 0.10204 # m
const crc3_rotor_radii = [10.20, 20.41, 30.61, 40.82, 51.02, 61.22, 66.33, 71.43, 76.53, 81.63, 86.73, 91.84, 96.94, 100.00] .* crc3_rotor_R / 100.0
const crc3_rotor_chord = [0.1195, 0.1340, 0.1607, 0.1920, 0.2142, 0.2200, 0.2190, 0.2100, 0.2039, 0.1911, 0.1758, 0.1581, 0.1272, 0.0982] .* crc3_rotor_R
const crc3_rotor_twist = [30.35612, 34.12497, 28.91561, 21.20027, 15.47515, 12.6852, 11.21356, 10.46311, 10.0537, 9.066763, 7.5527, 7.540607, 8.273884, 8.540983] .* pi/180

const crc3_rotor_x = 65.1925e-3 - 0.28717 * 0.3048/4 # m
const crc3_rotor_y = 242.55e-3/2 # m
const crc3_rotor_rhub = crc3_rotor_radii[1] # m
const crc3_rotor_rtip = crc3_rotor_radii[end] # m
const crc3_rotor_D = 2 * crc3_rotor_R
const crc3_rotor_cant = 10*pi/180 # inboard inclination of the rotor axes

const crc3_rotor_contours = [
    "e212.dat",
    "e212.dat",
    "e212.dat",
    "e212.dat",
    "e212.dat",
    "e212.dat",
    "e212.dat",
    "e212.dat",
    "e212.dat",
    "e212.dat",
    "e212.dat",
    "e212.dat",
    "e212.dat",
    "e212.dat"
]
# const crc3_rotor_contours = [
#     "contour10.20.csv",
#     "contour20.41.csv",
#     "contour30.61.csv",
#     "contour40.82.csv",
#     "contour51.02.csv",
#     "contour61.22.csv",
#     "contour66.33.csv",
#     "contour71.43.csv",
#     "contour76.53.csv",
#     "contour81.63.csv",
#     "contour86.73.csv",
#     "contour91.84.csv",
#     "contour96.94.csv",
#     "contour100.00.csv"
# ]
const crc3_rotor_polars = [
    "contour10.20_polar_5000.0_0.03.csv",
    "contour20.41_polar_5000.0_0.14.csv",
    "contour51.02_polar_5000.0_0.46.csv",
    "contour51.02_polar_5000.0_0.46.csv",
    "contour51.02_polar_5000.0_0.46.csv",
    "contour61.22_polar_5000.0_0.57.csv",
    "contour66.33_polar_5000.0_0.62.csv",
    "contour66.33_polar_5000.0_0.62.csv",
    "contour66.33_polar_5000.0_0.62.csv",
    "contour66.33_polar_5000.0_0.62.csv",
    "contour66.33_polar_5000.0_0.62.csv",
    "contour66.33_polar_5000.0_0.62.csv",
    "contour96.94_polar_5000.0_0.95.csv",
    "contour100.00_polar_5000.0_1.0.csv"
]

const crc3_wing_spacing = 0.254 # distance from lower wing to upper wing
const crc3_wing_b = 20 * 0.0254 # m
const crc3_wing_c = 0.28717 * 0.3048 # m

const crc3_mass = 3 / 2.2 # kg
const crc3_JInertia = 0.0111111 # kg-m^2
