from math import sqrt
import numpy as np
import simulation

# Functions to handle the geometry of a cylindrical elbow pipe with radius
# 1 and dimensions A, B
#
#             B
#       /----------/
#     +-------------     +--> y
#   \ | O                |
#   | |   +---------     v
# A | |   |              x
#   | |   |
#   \ |   |
#
# The junction of the pipe lies at the origin O.

# A "large distance" used to indicate there's no collision
Sm = 1e99
eps = 1e-8


# Generates a new particle at opening A with a random direction
def newParticle(A):
    while True:
        yz = 2 * np.random.rand(2) - 1
        if yz.dot(yz) < 1:
            break
    return (np.array([A, yz[0], yz[1]]),
            simulation.newDirection(openingABasis))


openingABasis = np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]])


# Returns a local Cartesian basis of the surface of pipe A such that the third
# axis is the normal
def pipeABasis(pos):
    return np.array([[1, 0, 0], [0, -pos[2], -pos[1]],
                     [0, pos[1], -pos[2]]]) / np.linalg.norm(pos[1:3])


permXY = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]])


# Pipe B is identical if we just permute the x and y components
def pipeBBasis(pos):
    return permXY.dot(pipeABasis(permXY.dot(pos)))


# Functions to calculate distances to the surfaces
# p refers to the position vector and s direction vector

# Distance to intersection with the cylinder y^2 + z^2 = 1
def distanceToPipeA(p, s):
    s12, p12 = s[1:3], p[1:3]
    a, b, c = s12.dot(s12), s12.dot(p12), p12.dot(p12) - 1
    bac = b**2 - a * c
    if abs(a) < eps or bac < 0:
        return Sm
    S = (sqrt(bac) - b) / a
    pNew = p + S * s
    # Check if we're on the right side of the plane x = y
    if pNew[0] < pNew[1]:
        return Sm
    return S


# Distance to intersection with the plane x = A
def distanceToOpeningA(p, s, A):
    if abs(s[0]) < eps:
        return Sm
    return (A - p[0]) / s[0]


# If the distance is negative (backwards) or very small, exclude it by setting
# it to a large value
def handleDistance(S):
    if S < eps:
        return Sm
    return S


# Computes the distance to the next collision. Returns a tuple with the
# distance and index of the collision:
# -1: no collision
# 0: opening A
# 1: opening B
# 2: pipe A
# 3: pipe B
def nextCollision(p, s, A, B):
    pB, sB = permXY.dot(p), permXY.dot(s)
    distances = [handleDistance(S) for S in [
        distanceToOpeningA(p, s, A), distanceToOpeningA(pB, sB, B),
        distanceToPipeA(p, s), distanceToPipeA(pB, sB)
    ]]
    Smin = min(distances)
    idx = distances.index(Smin)
    if Smin == Sm:
        idx = -1

    pNew = p + Smin * s
    if idx == 2:
        sNew = simulation.newDirection(pipeABasis(pNew))
    elif idx == 3:
        sNew = simulation.newDirection(pipeBBasis(pNew))
    else:
        sNew = None

    return Smin, idx, pNew, sNew
