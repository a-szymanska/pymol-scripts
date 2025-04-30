# PyMOL Script Collection  
<!-- ![License](https://img.shields.io/badge/license-MIT-green)  -->

This repository is a small collection of PyMOL scripts created for visualization and analysis.

## Overview  
`mark_termini.py`  
Marks all C-termini or N-termini in a protein
```python
mark_termini(name="6RVV", term="C")
```

\
`pack_cgo.py`  
Fits a CGO solid to a molecular structure, possible bounding shapes: sphere, cylinder, cuboid, cone
```python
pack_sphere("1HZH", color=(1,0,1), margin=10)
pack_cuboid("1HZH", cube=False)
```

\
`align_structures.py`  
Aligns structure A at a specified point and orients it perpendicular to structure B based on provided points
```python
align(v0, v1, p0, p1, p2, "1HZH")
```
Points v0, v1 determine the axis along structure A, points p0, p1, p2 determine the perpendicular plane on structure B. Structure A is 'attached' to p0 at v0.

\
`measure.py`  
Generates a bounding box for a structure and measures its size in the three dimensions
```python
m = Measurement("1HZH") # the default is "all"
m.rotate('x', -10)
print(m)
```


\
`min_distance.py`  
Finds minimum distance between atoms of the two input structures. It uses a Divide-and-Conquer algorithm with a complexity of $O(n \log n)$, where $n$ is the total number of atoms, and thus it runs several dozen times faster than the brute force approach.
```python
dist = min_dist("1AAC", "5AAC")
print(dist)
```

## Dependencies
- PyMOL
- NumPy
- SciPy
