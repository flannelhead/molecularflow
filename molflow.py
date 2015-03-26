import sys
import numpy as np
import elbow
import simulation

if len(sys.argv) < 5:
    print('Usage: python molflow.py A B N dx <outfile>')
    sys.exit(0)

# Dimensions of the elbow tube
A = float(sys.argv[1])
B = float(sys.argv[2])

# Number of tracked particles
N = int(sys.argv[3])
# Maximum number of collisions before we stop tracking a particle
MAX_COLLISIONS = 30

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
    p, s = elbow.newParticle(A)
    remainder = 0

    for j in range(MAX_COLLISIONS):
        distance, idx, pNew, sNew = elbow.nextCollision(p, s, A, B)

        # Stop if there's no next collision
        if idx == -1:
            nRejected += 1
            break

        remainder = simulation.traceSegment(distance, remainder, dx, p, s, nArr,
            origin, nPoints)

        # Break out if the particle has exited
        if idx == 0:
            break
        elif idx == 1:
            nB += 1
            break

        p, s = pNew, sNew

print('{0} particles transmitted, {1} particles rejected'.format(nB,
    nRejected))
print('Transmission probability Pr = ' + str(nB / (N - nRejected)))

# Save computation results to a file
if len(sys.argv) > 5: filename = sys.argv[5]
else: filename = 'output.npz'
f = open(filename, 'wb')
np.savez(f, X = grid[0], Y = grid[1], Z = grid[2], C = nArr)
f.close()
print('Heatmap saved to ' + filename + '.')

