from pymol import cgo
import numpy as np
from scipy.spatial.distance import euclidean

def pack_sphere(name : str):
    coord_center = cmd.centerofmass(name)
    model = cmd.get_model(name)
    coords = [atom.coord for atom in model.atom]
    distances = [euclidean(coord_center, coord) for coord in coords]
    x, y, z, r = *coord_center, np.max(distances) + 1
    sphere = [cgo.SPHERE, x, y, z, r]
    cmd.load_cgo(sphere, f"{name}_sphere")
    cmd.set("cgo_transparency", 0.6, f"{name}_sphere")


def pack_cylinder(name : str):
    cmd.orient()
    model = cmd.get_model(name)
    coords = np.array([atom.coord for atom in model.atom])
    coord_range = [np.ptp(coords[:, i]) for i in range(3)]
    idx = np.argmax(coord_range)
    pos0 = [0, 0, 0]
    pos1 = [0, 0, 0]
    pos0[idx] = np.min(coords[:, idx])
    pos1[idx] = np.max(coords[:, idx])
    r = 0
    for i in range(3):
        if i == idx:
            continue
        r = max(r, coord_range[i]/2)
        pos0[i] = pos1[i] = (np.min(coords[:, i]) + np.max(coords[:, i])) / 2
    r += 5
    white = [1, 1, 1]
    cylinder = [cgo.CONE, *pos0, *pos1, r, r, *white, *white, 0, 1]
    cmd.load_cgo(cylinder, f"{name}_cylinder")
    cmd.set("cgo_transparency", 0.6, f"{name}_cylinder")


def pack_cuboid(name : str, cube : bool = False):
    pass


def pack_cone(name : str, truncated : bool = False):
    pass