from pymol import cmd, cgo
import numpy as np

class Measurement:
    def __init__(self, name: str = "all"):
        self.name = name
        self.coords = []
        self.fit()
        self.measure()

    def fit(self):
        model = cmd.get_model(self.name)
        points = np.array([atom.coord for atom in model.atom])
        min_coords = np.min(points, axis=0)
        max_coords = np.max(points, axis=0)
        self.coords = [min_coords, max_coords]

    def measure(self):
        x_min, y_min, z_min = self.coords[0]
        x_max, y_max, z_max = self.coords[1]
        cuboid = [
            cgo.BEGIN, cgo.LINES,  
            cgo.VERTEX, x_min, y_min, z_max, cgo.VERTEX, x_max, y_min, z_max,
            cgo.VERTEX, x_min, y_min, z_max, cgo.VERTEX, x_min, y_max, z_max,
            cgo.VERTEX, x_min, y_min, z_max, cgo.VERTEX, x_min, y_min, z_min,
            cgo.END
        ]
        cmd.delete(f"{self.name}_measurement")
        cmd.delete(f"{self.name}_labels")
        cmd.load_cgo(cuboid, f"{self.name}_measurement")
        
        labels = [
            ((x_min + x_max) / 2, y_min, z_max, f"{abs(x_max - x_min):.2f}"),
            (x_min, (y_min + y_max) / 2, z_max, f"{abs(y_max - y_min):.2f}"),
            (x_min, y_min, (z_min + z_max) / 2, f"{abs(z_max - z_min):.2f}")
        ]
        for x, y, z, text in labels:
            cmd.pseudoatom(f"{self.name}_labels", pos=[x, y, z], label=text)
        cmd.zoom()

    def rotate(self, axis: str, angle: int):
        if axis.lower() not in ['x', 'y', 'z']:
            print('Rotation only possible in axis x, y or z.')
            return
        cmd.rotate(axis, angle, self.name)
        self.fit()
        self.measure()

    def __repr__(self):
        x_min, y_min, z_min = self.coords[0]
        x_max, y_max, z_max = self.coords[1]
        return (f"({x_max - x_min:.2f}, {y_max - y_min:.2f}, {z_max - z_min:.2f})")
