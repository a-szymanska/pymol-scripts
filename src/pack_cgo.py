from pymol import cmd, cgo
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean
from scipy.optimize import minimize

TRANSPARENCY = 0.4
SPHERE_QUALITY = 4
CONE_QUALITY = 100
LINE_WIDTH = 4.0
COLOR = (1, 1, 1)
EPSILON = 1e-6
MAX_ITER = 5


# ----------------------- Helper functions -----------------------
def get_available_name(base_name):
    suffix = ''
    existing_objects = cmd.get_names("all")
    while f"{base_name}{suffix}" in existing_objects:
        suffix = (suffix or 0) + 1
    return f"{base_name}{suffix}" if suffix else base_name


# p - point on the axis corresponding to x, v - axis vector
def get_distance_to_axis(x, p, v):
    return np.linalg.norm(np.cross(v, x - p)) / np.linalg.norm(v)


def find_base_centers(points, p, v):
    projections = [np.dot((x - p), v) for x in points]
    p_min = min(projections)
    p_max = max(projections)
    p1 = p + p_min * v
    p2 = p + p_max * v
    return p1, p2


# p - point on the axis, v - axis vector, r - cylinder radius
def f_objective_cylinder(params, points):
    p, v, r = params[:3], params[3:6], params[6]
    v /= np.linalg.norm(v)
    cost = 0
    for x in points:
        dist = get_distance_to_axis(x, p, v)
        cost += (dist - r)**2
    return cost


def fit_cylinder(points):
    centroid = np.mean(points, axis=0)
    v_init = np.array([1.0, 0.0, 0.0])
    r_init = np.mean(np.linalg.norm(points - centroid, axis=1)) / 2
    params = np.concatenate([centroid, v_init, [r_init]])
    result = minimize(f_objective_cylinder, params, args=(points,), method='L-BFGS-B', options={'maxiter': MAX_ITER})
    p_min = result.x[:3]
    v_min = result.x[3:6]
    r_min = result.x[6]
    return p_min, v_min / np.linalg.norm(v_min), r_min


def get_line_cuboid(min_coords, max_coords):
    x_min, y_min, z_min = min_coords
    x_max, y_max, z_max = max_coords
    cuboid = [
        cgo.BEGIN, cgo.LINES, cgo.COLOR, *COLOR,  

        cgo.VERTEX, x_min, y_min, z_min, cgo.VERTEX, x_max, y_min, z_min,
        cgo.VERTEX, x_max, y_min, z_min, cgo.VERTEX, x_max, y_max, z_min,
        cgo.VERTEX, x_max, y_max, z_min, cgo.VERTEX, x_min, y_max, z_min,
        cgo.VERTEX, x_min, y_max, z_min, cgo.VERTEX, x_min, y_min, z_min,

        cgo.VERTEX, x_min, y_min, z_max, cgo.VERTEX, x_max, y_min, z_max,
        cgo.VERTEX, x_max, y_min, z_max, cgo.VERTEX, x_max, y_max, z_max,
        cgo.VERTEX, x_max, y_max, z_max, cgo.VERTEX, x_min, y_max, z_max,
        cgo.VERTEX, x_min, y_max, z_max, cgo.VERTEX, x_min, y_min, z_max,

        cgo.VERTEX, x_min, y_min, z_min, cgo.VERTEX, x_min, y_min, z_max,
        cgo.VERTEX, x_max, y_min, z_min, cgo.VERTEX, x_max, y_min, z_max,
        cgo.VERTEX, x_max, y_max, z_min, cgo.VERTEX, x_max, y_max, z_max,
        cgo.VERTEX, x_min, y_max, z_min, cgo.VERTEX, x_min, y_max, z_max,

        cgo.END
    ]
    return cuboid


def get_colored_cuboid(min_coords, max_coords, color):
    x_min, y_min, z_min = min_coords
    x_max, y_max, z_max = max_coords

    cuboid = [
        cgo.BEGIN, cgo.TRIANGLES,
        cgo.COLOR, *color,

        cgo.VERTEX, x_min, y_min, z_min,
        cgo.VERTEX, x_max, y_min, z_min,
        cgo.VERTEX, x_max, y_max, z_min,
        cgo.VERTEX, x_max, y_max, z_min,
        cgo.VERTEX, x_min, y_max, z_min,
        cgo.VERTEX, x_min, y_min, z_min,

        cgo.VERTEX, x_min, y_min, z_max,
        cgo.VERTEX, x_max, y_min, z_max,
        cgo.VERTEX, x_max, y_max, z_max,
        cgo.VERTEX, x_max, y_max, z_max,
        cgo.VERTEX, x_min, y_max, z_max,
        cgo.VERTEX, x_min, y_min, z_max,

        cgo.VERTEX, x_min, y_min, z_min,
        cgo.VERTEX, x_max, y_min, z_min,
        cgo.VERTEX, x_max, y_min, z_max,
        cgo.VERTEX, x_max, y_min, z_max,
        cgo.VERTEX, x_min, y_min, z_max,
        cgo.VERTEX, x_min, y_min, z_min,

        cgo.VERTEX, x_min, y_max, z_min,
        cgo.VERTEX, x_max, y_max, z_min,
        cgo.VERTEX, x_max, y_max, z_max,
        cgo.VERTEX, x_max, y_max, z_max,
        cgo.VERTEX, x_min, y_max, z_max,
        cgo.VERTEX, x_min, y_max, z_min,

        cgo.VERTEX, x_min, y_min, z_min,
        cgo.VERTEX, x_min, y_max, z_min,
        cgo.VERTEX, x_min, y_max, z_max,
        cgo.VERTEX, x_min, y_max, z_max,
        cgo.VERTEX, x_min, y_min, z_max,
        cgo.VERTEX, x_min, y_min, z_min,

        cgo.VERTEX, x_max, y_min, z_min,
        cgo.VERTEX, x_max, y_max, z_min,
        cgo.VERTEX, x_max, y_max, z_max,
        cgo.VERTEX, x_max, y_max, z_max,
        cgo.VERTEX, x_max, y_min, z_max,
        cgo.VERTEX, x_max, y_min, z_min,

        cgo.END
    ]
    return cuboid


# ------------------------ Main functions ------------------------
def pack_sphere(name : str, color : tuple = COLOR, margin : float = 1.0):
    if margin < 0:
        print("! Warning: negative margin !")

    coord_center = cmd.centerofmass(name)
    model = cmd.get_model(name)
    coords = [atom.coord for atom in model.atom]
    distances = [euclidean(coord_center, coord) for coord in coords]
    x, y, z, r = *coord_center, np.max(distances) + margin
    sphere = [cgo.COLOR, *color, cgo.SPHERE, x, y, z, r]

    name_packed = get_available_name(f"{name}_sphere")
    print(name_packed)
    cmd.load_cgo(sphere, name_packed)
    cmd.set("cgo_transparency", TRANSPARENCY, name_packed)
    cmd.set("cgo_sphere_quality", SPHERE_QUALITY)


def pack_cylinder(name : str, color : tuple = COLOR, margin : float = 1.0):
    if margin < 0:
        print("! Warning: negative margin !")

    model = cmd.get_model(name)
    points = np.array([atom.coord for atom in model.atom])

    p0_fitted, axis_vector, r = fit_cylinder(points)

    distances = [get_distance_to_axis(x, p0_fitted, axis_vector) for x in points]
    max_dist = np.max(distances)
    r = max(r, max_dist) + EPSILON + margin
    p1_base, p2_base = find_base_centers(points, p0_fitted, axis_vector)
    axis_vector = axis_vector / np.linalg.norm(axis_vector)  # Normalize axis vector
    p1_base -= axis_vector * margin
    p2_base += axis_vector * margin
    
    cylinder = [cgo.CONE, *p1_base, *p2_base, r, r, *color, *color, 0, 1]

    name_packed = get_available_name(f"{name}_cylinder")
    cmd.load_cgo(cylinder, name_packed)
    cmd.set("cgo_transparency", TRANSPARENCY, name_packed)
    cmd.set("cone_quality", CONE_QUALITY)


def pack_cuboid(name: str, cube: bool = False, color : tuple = None, margin: float = 1.0):
    if margin < 0:
        print("! Warning: negative margin !")

    model = cmd.get_model(name)
    points = np.array([atom.coord for atom in model.atom])

    min_coords = np.min(points, axis=0) - margin
    max_coords = np.max(points, axis=0) + margin

    if cube:
        max_range = np.max(max_coords - min_coords)
        center = (min_coords + max_coords) / 2
        min_coords = center - max_range / 2
        max_coords = center + max_range / 2
    if color is None:
        cuboid = get_line_cuboid(min_coords, max_coords)
    else:
        cuboid = get_colored_cuboid(min_coords, max_coords, color)

    name_packed = get_available_name(f"{name}_{'cube' if cube else 'cuboid'}")
    cmd.load_cgo(cuboid, name_packed)
    cmd.set("cgo_transparency", TRANSPARENCY, name_packed)
