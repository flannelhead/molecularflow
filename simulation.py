import random as rnd
import numpy as np
from math import pi, sqrt, sin, cos, asin
from scipy.stats import maxwell

a = 8 / pi


# Generate a random scattering angle
def cosineDirection():
    phi = 2.0 * pi * rnd.random()
    theta = asin(sqrt(rnd.random()))
    return np.array([sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)])


# Generate a speed from the Maxwell-Boltzmann distribution
def mbSpeed(T, kOverM, sampleMB=False):
    if sampleMB:
        return maxwell.ppf(rnd.random(), scale=sqrt(kOverM*T))
    else:
        return sqrt(a*kOverM*T)


# Generate a scattering direction with respect to a given surface
def newDirection(basis):
    return basis.dot(cosineDirection())


# Create a cube grid for recording the particle locations
def cubeGrid(dimensions, origin, cellSize):
    nPoints = (dimensions / cellSize).astype(int)
    arr = np.zeros(nPoints)
    grid = [np.linspace(0, dimensions[i], nPoints[i] + 1) + origin[i]
            for i in [0, 1, 2]]
    return arr, grid, nPoints


# Trace a segment of a particle's travel
def traceSegment(distance, remainder, dx, dt, p, s, v, arr, offset, nPoints):
    sRemainder = v * remainder
    if distance - sRemainder < 0:
        return remainder - distance / v
    ds = v * dt
    nSteps = int((distance - sRemainder) / ds)
    sampleDistances = np.linspace(0, nSteps * ds, nSteps + 1) + sRemainder
    for S in sampleDistances:
        i = np.maximum(np.minimum(np.floor((p + S * s - offset) / dx),
                       nPoints - 1), 0)
        arr[i[0], i[1], i[2]] += 1
    return ((nSteps + 1) * ds + sRemainder - distance) / v
