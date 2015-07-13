import random as rnd
import numpy as np
from math import pi, sqrt, sin, cos, asin
from scipy.stats import maxwell


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
        return maxwell.mean(scale=sqrt(kOverM*T))


# Generate a scattering direction with respect to a given surface
def newDirection(basis):
    return basis.dot(cosineDirection())


# Create a cube grid for recording the particle locations
def cubeGrid(dimensions, origin, cellSize, steps, dt):
    nPoints = (dimensions / cellSize).astype(int)
    grid = [np.linspace(0, dimensions[i], nPoints[i] + 1) + origin[i]
            for i in [0, 1, 2]]
    steps = max(steps, 1)
    arr = np.zeros(np.insert(nPoints, 0, steps))
    return arr, grid, nPoints


# Trace a segment of a particle's travel
def traceSegment(distance, remainder, gridSize, dt, p, s, v, arr, offset,
                 nPoints, steps, currentStep, isSteadyState):
    ds = v * dt
    S = v * remainder
    while S < distance:
        i = np.maximum(np.minimum(np.floor((p + S * s - offset) / gridSize),
                       nPoints - 1), 0)
        if isSteadyState:
            arr[0, i[0], i[1], i[2]] += 1
        else:
            if currentStep >= steps:
                break
            arr[currentStep, i[0], i[1], i[2]] += 1
        currentStep += 1
        S += ds
    return (S - distance) / v, currentStep
