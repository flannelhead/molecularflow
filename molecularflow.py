import sys
import numpy as np
from scipy import constants
import elbow
import simulation

if len(sys.argv) < 5:
    print('Usage: python molecularflow.py A B N dx <outfile>')
    sys.exit(0)

T = 293  # room temperature
km = constants.k / (28 * constants.value('atomic mass constant'))  # nitrogen
Q = 0.01  # Flow rate

# Dimensions of the elbow tube
A = float(sys.argv[1])
B = float(sys.argv[2])

# Number of tracked particles
N = int(sys.argv[3])
# Maximum number of collisions before we stop tracking a particle
MAX_COLLISIONS = 50
dt = 0.0003

# Define the problem domain
origin = np.array([-1, -1, -1])
dimensions = np.array([A + 1, B + 1, 2])
dx = float(sys.argv[4])
nArr, grid, nPoints = simulation.cubeGrid(dimensions, origin, dx)

# Number of transmitted particles
nB = 0
nRejected = 0

print('Computing with A={0}, B={1}, N={2}'.format(A, B, N))

for i in range(N):
    p, s, v = elbow.newParticle(A, T, km)
    remainder = 0

    for j in range(MAX_COLLISIONS):
        distance, idx, pNew, sNew, vNew = elbow.nextCollision(p, s,
                                                              A, B, T, km)

        # Stop if there's no next collision
        if idx == -1:
            nRejected += 1
            break

        remainder = simulation.traceSegment(distance, remainder, dx, dt, p,
                                            s, v, nArr, origin, nPoints)

        # Break out if the particle has exited
        if idx == 0:
            break
        elif idx == 1:
            nB += 1
            break

        p, s, v = pNew, sNew, vNew

print('{0} particles transmitted, {1} particles rejected'
      .format(nB, nRejected))
print('Transmission probability Pr = ' + str(nB / (N - nRejected)))

nArr *= Q*dt / (N * dx**3 * constants.k * T)

# Save the computation results to a file
if len(sys.argv) > 5:
    filename = sys.argv[5]
else:
    filename = 'output.npz'
f = open(filename, 'wb')
np.savez(f, X=grid[0], Y=grid[1], Z=grid[2], C=nArr)
f.close()
print('Heatmap saved to ' + filename + '.')
