import numpy as np
from scipy import constants
from . import elbow
from . import simulation


def runSimulation(A, B, Q=0.01, Z=28, T=300, N=10000, dx=0.2, dt=1e-3,
                  MAX_COLLISIONS=50, sampleMB=False, filename=None):
    km = constants.k / (Z * constants.value('atomic mass constant'))

    # Define the problem domain
    offset = np.array([-1, -1, -1])
    dimensions = np.array([A + 1, B + 1, 2])
    nArr, grid, nPoints = simulation.cubeGrid(dimensions, offset, dx)

    # Number of transmitted particles
    nB = 0
    nRejected = 0

    print('Computing with A={0}, B={1}, N={2}'.format(A, B, N))

    for i in range(N):
        p, s = elbow.newParticle(A)
        v = simulation.mbSpeed(T, km, sampleMB)
        remainder = 0

        for j in range(MAX_COLLISIONS):
            distance, idx, pNew, sNew = elbow.nextCollision(p, s, A, B)
            vNew = simulation.mbSpeed(T, km, sampleMB)

            # Stop if the particle was rejected
            if idx == -1:
                nRejected += 1
                break

            remainder = simulation.traceSegment(distance, remainder, dx, dt, p,
                                                s, v, nArr, offset, nPoints)

            # Break out if the particle has exited
            if idx == 0:
                break
            elif idx == 1:
                nB += 1
                break

            p, s, v = pNew, sNew, vNew

    nArr *= Q*dt / (N * dx**3 * T * constants.k)
    Pr = nB / (N - nRejected)

    print('{0} particles transmitted, {1} particles rejected'
          .format(nB, nRejected))
    print('Transmission probability Pr = ' + str(Pr))

    # Save the computation results to a file
    if filename is not None:
        np.savez(filename, X=grid[0], Y=grid[1], Z=grid[2], C=nArr,
                 Pr=np.array([Pr]))
        print('Results saved to ' + filename + '.')

    return Pr, grid, nArr
