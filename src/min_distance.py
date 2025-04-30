from pymol import cmd, cgo
from scipy.spatial.distance import euclidean as dist

def strip_closest(strip, d):
    min_d = d
    min_points = None
    strip.sort(key=lambda point: point[0][1])
    for i in range(len(strip)):
        for j in range(i + 1, len(strip)):
            if strip[i][1] != strip[j][1] and dist(strip[i][0], strip[j][0]) < min_d:
                min_d = dist(strip[i][0], strip[j][0])
                min_points = (strip[i][0], strip[j][0])
            else:
                break
    return min_d, min_points


def min_dist_helper(points, l, r):
    if r - l <= 10:
        min_d = float('inf')
        min_points = None
        for i, p in enumerate(points[l: r]):
            for q in points[i + 1: r]:
                if p[1] != q[1] and dist(p[0], q[0]) < min_d:
                    min_d = dist(p[0], q[0])
                    min_points = (p[0], q[0])
        return min_d, min_points

    mid = (l + r) // 2
    mid_x = points[mid][0][0]
    d_left, points_left = min_dist_helper(points, l, mid)
    d_right, points_right = min_dist_helper(points, mid, r)
    if d_left < d_right:
        d = d_left
        min_points = points_left
    else:
        d = d_right
        min_points = points_right
    strip = [points[i] for i in range(l, r) if abs(points[i][0][0] - mid_x) < d]
    d_center, points_center = strip_closest(strip, d)
    if (d <= d_center):
        return d, min_points
    else:
        return d_center, points_center


def min_dist(objA : str, objB : str) -> float:
    pointsA = cmd.get_coords(objA, 1)
    pointsB = cmd.get_coords(objB, 1)
    pointsAB = [(p, 0) for p in pointsA] + [(p, 1) for p in pointsB]
    pointsAB.sort(key=lambda p: p[0][0])
    min_d, (pointA, pointB) = min_dist_helper(pointsAB, 0, len(pointsAB))
    for i, point in enumerate([pointA, pointB]):
        print(point)
        sphere = [cgo.SPHERE, *point, 0.4]
        cmd.load_cgo(sphere, f"point{i}")
        cmd.pseudoatom("tmp", pos=list(point))
        cmd.select("nearest", "all within 5 of tmp")
        cmd.show("licorice", "nearest")
        cmd.delete("tmp")
        cmd.delete("nearest")
    line = [
        cgo.BEGIN, cgo.LINES,
        cgo.VERTEX, *pointA,
        cgo.VERTEX, *pointB,
        cgo.END
    ]
    cmd.load_cgo(line, "dist_line")
    mid = [(a + b) / 2 for a, b in zip(pointA, pointB)]
    label_pos = [mid[0], mid[1], mid[2] + 1.0]
    cmd.pseudoatom(object="dist_label", pos=label_pos, label=f"{min_d:.2f}")
    cmd.zoom()
    print(f"Minimum distance: {min_d:.4f}")
    return min_d