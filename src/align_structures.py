from pymol import cmd
import numpy as np

'''
v0, v1, p0, p1, p2 - selections, may consist of several atoms
name - name of structure A (the one to be aligned)

v0, v1 determine the axis along structure A
p0, p1, p2 determine the perpendicular plane on structure B
v0 is 'attached' to p0
'''
def align(v0, v1, p0, p1, p2, name):

    def get_coord(selection):
        model = cmd.get_model(selection)
        return np.array(model.atom[0].coord)

    # Coordinates approximated by choosing
    # one atom from each selection
    v0_coord = get_coord(v0)
    v1_coord = get_coord(v1)
    p0_coord = get_coord(p0)
    p1_coord = get_coord(p1)
    p2_coord = get_coord(p2)

    v_axis = v1_coord - v0_coord
    p_x = p1_coord - p0_coord
    p_y = p2_coord - p0_coord

    p_normal = np.cross(p_x, p_y)
    p_normal /= np.linalg.norm(p_normal)
    v_axis = v_axis - np.dot(v_axis, p_normal) * p_normal
    v_axis /= np.linalg.norm(v_axis)
    v_normal = np.cross(v_axis, p_normal)
    v_normal /= np.linalg.norm(v_normal)

    # Transformation matrices
    A = np.column_stack([v_axis, v_normal, np.cross(v_axis, v_normal)])
    B = np.column_stack([p_x / np.linalg.norm(p_x), p_y / np.linalg.norm(p_y), p_normal])
    R = B @ np.linalg.inv(A)

    cmd.alter_state(1, name, "(x,y,z) = np.dot(R, np.array([x,y,z]) - v0_coord) + p0_coord", space={'np': np, 'R': R, 'v0_coord': v0_coord, 'p0_coord': p0_coord})
    cmd.refresh()
    print("If the structure orientation is the opposite of what is expected, call the function again with the same arguments. That should do the trick.")
