from pymol import cmd, cgo
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean
from scipy.optimize import curve_fit
from copy import copy

TRANSPARENCY = 0.4
SPHERE_QUALITY = 4
CONE_QUALITY = 100
COLOR = (1, 1, 1)
EPSILON = 1e-6


def pack_sphere(name : str, margin : float = 1.0):
    if margin < 0:
        print("! Warning: negative margin !")
    coord_center = cmd.centerofmass(name)
    model = cmd.get_model(name)
    coords = [atom.coord for atom in model.atom]
    distances = [euclidean(coord_center, coord) for coord in coords]
    x, y, z, r = *coord_center, np.max(distances) + margin
    sphere = [cgo.SPHERE, x, y, z, r]
    cmd.load_cgo(sphere, f"{name}_sphere")
    cmd.set("cgo_transparency", TRANSPARENCY, f"{name}_sphere")
    cmd.set("cgo_sphere_quality", SPHERE_QUALITY)


def pack_cylinder(name : str, margin : float = 1.0):
    if margin < 0:
        print("! Warning: negative margin !")
    cmd.orient(name)
    coord_center = cmd.centerofmass(name)
    model = cmd.get_model(name)
    coords = np.array([atom.coord for atom in model.atom])
    axis_idx = np.argmax(np.ptp(coords, axis=0))
    max_dist = 0
    for coord in coords:
        coord_center_axis = coord_center
        coord_center_axis[axis_idx] = coord[axis_idx]
        max_dist = max(max_dist, euclidean(coord, coord_center_axis))
    r = max_dist + margin
    pos0 = copy(coord_center)
    pos1 = copy(coord_center)
    pos0[axis_idx] = np.min(coords[:, axis_idx]) + margin
    pos1[axis_idx] = np.max(coords[:, axis_idx]) + margin
    cylinder = [cgo.CONE, *pos0, *pos1, r, r, *COLOR, *COLOR, 0, 1]
    cmd.load_cgo(cylinder, f"{name}_cylinder")
    cmd.set("cgo_transparency", TRANSPARENCY, f"{name}_cylinder")
    cmd.set("cone_quality", CONE_QUALITY)


def pack_cone(name : str, truncated : bool = False, margin : float = 1.0):
    cmd.orient(name)
    if margin < 0:
        print("! Warning: negative margin !")
    cmd.orient(name)
    coord_center = cmd.centerofmass(name)
    model = cmd.get_model(name)
    coords = np.array([atom.coord for atom in model.atom])
    axis_idx = np.argmax(np.ptp(coords, axis=0))
    distances = []
    for coord in coords:
        coord_center_axis = coord_center
        coord_center_axis[axis_idx] = coord[axis_idx]
        distances.append(euclidean(coord, coord_center_axis))
    x_data = coords[:, axis_idx]
    y_data = distances
    popt, pcov = curve_fit((lambda x,a,b,c: a*x+b), x_data, y_data)
    a, b = popt[:2]
    f = lambda x: a*x + b
    shift = np.max(y_data - f(x_data)) + EPSILON
    b += shift
    r_big = f(np.max(coords[:, axis_idx])) + margin
    pos_base = copy(coord_center)
    pos_base[axis_idx] = np.max(coords[:, axis_idx]) + margin
    pos_top = copy(coord_center)
    if truncated:
        r_small = f(np.min(coords[:, axis_idx])) + margin
        pos_top[axis_idx] = np.min(coords[:, axis_idx]) + margin
    else:
        r_small = EPSILON
        pos_top[axis_idx] = -b/a
    cone = [cgo.CONE, *pos_top, *pos_base, r_small, r_big, *COLOR, *COLOR, 0, 1]
    cmd.load_cgo(cone, f"{name}_cone")
    cmd.set("cgo_transparency", TRANSPARENCY, f"{name}_cone")
    cmd.set("cone_quality", CONE_QUALITY)


def pack_cuboid(name : str, cube : bool = False, margin : float = 1.0):
    if margin < 0:
        print("! Warning: negative margin !")
    pass
