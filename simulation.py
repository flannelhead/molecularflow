import random as rnd
import numpy as np
from math import pi, sqrt, sin, cos, asin

# Generate a random scattering angle
def lambertDirection():
    phi = 2.0 * pi * rnd.random()
    theta = asin(sqrt(rnd.random()))
    return np.array([sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)])

# Generate a scattering direction with respect to a given surface
def newDirection(basis):
    return basis.dot(lambertDirection())

# Create a cube grid for recording the particle locations
def cubeGrid(dimensions, origin, cellSize):
    nPoints = (dimensions / cellSize).astype(int)
    arr = np.zeros(nPoints)
    grid = [np.linspace(0, dimensions[i], nPoints[i] + 1) + origin[i]
            for i in [0, 1, 2]]
    return arr, grid, nPoints

# Trace a segment of a particle's travel
def traceSegment(distance, remainder, dx, p, s, arr, origin, nPoints):
    nSteps = int((distance - remainder) / dx)
    sampleDistances = np.linspace(0, nSteps * dx, nSteps + 1) + remainder
    remainder = (nSteps + 1) * dx + remainder - distance
    for S in sampleDistances:
        i = np.maximum(np.minimum(np.floor((p + S * s - origin) / dx),
            nPoints - 1), 0)
        arr[i[0], i[1], i[2]] += 1
    return remainder

