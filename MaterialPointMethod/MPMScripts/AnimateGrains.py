import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import sys

if len(sys.argv) < 4:
    print("Usage: AnimateGrains.py <file name> <x dim> <y dim>")
    exit()

outfile = open(sys.argv[1])
x_dim = int(sys.argv[2])
y_dim = int(sys.argv[3])

num_steps  = int(outfile.readline())

num_points = int(outfile.readline())

time_points = np.zeros((num_steps, num_points, 2))

timestep = -1
point = -1
axis = 0

for line in outfile:
    line = line.strip()
    if "Timestep" in line:
        timestep += 1
        point = -1
    elif "Point" in line:
        point += 1
    elif len(line) > 0:
        time_points[timestep][point][axis] = float(line)
        axis = (axis + 1) % 2

fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(7, 6.5)

ax = plt.axes(xlim=(-x_dim, x_dim), ylim=(-y_dim, y_dim))
patches = []
for p in range(num_points):
    patches.append(plt.Circle((time_points[0][p][0], 
                            time_points[0][p][1]), 0.25, fc='y'))

def init():
    for p in range(num_points):
        patches[p].center = (time_points[0][p][0], 
                            time_points[0][p][1])
        ax.add_patch(patches[p])
    return patches,

def animate(i):
    for p in range(num_points):
        patches[p].center = (time_points[i][p][0], 
                            time_points[i][p][1])
    return patches,

anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=num_steps, 
                               interval=1,
                               blit=False)

plt.show()