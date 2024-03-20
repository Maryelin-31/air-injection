# -*- coding: utf-8 -*-
from zmlx import *

# Create seepage mesh.
def create_mesh():
    mesh = SeepageMesh.create_cube(np.linspace(0, 30, 121), np.linspace(0, 15, 61), (-1.5, 1.5))
    # mesh = SeepageMesh.create_cube(np.linspace(0, 50, 101), np.linspace(0, 50, 26), np.linspace(0, 10, 10))
    return mesh
