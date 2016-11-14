import PIL
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib._png import read_png
from scipy.interpolate import interp1d

# Setup 3D Plot --------------------------------------------------------
fig = plt.figure(figsize=(10, 10)) 
ax = fig.add_subplot(111, projection='3d')

ax.set_xlim([-2, 2])
ax.set_ylim([-2, 2])
ax.set_zlim([-2, 2])

ax.set_xlabel('X Dim')
ax.set_ylabel('Y Dim')
ax.set_zlabel('Z Dim')

# Plot Earth -----------------------------------------------------------
# The approach is to load the blue marble image into the an image of RGB
# tuples, with each of R, G, B being in [0, 1]. This is used as the
# coloring for a surface plot of a sphere.
bm = PIL.Image.open('bluemarble.jpg')
bm = np.array(bm.resize([d/5 for d in bm.size]))
bm = bm.astype(float) / 257.

lons = np.linspace(-180, 180, bm.shape[1])       * np.pi / 180 
lats = np.linspace( -90,  90, bm.shape[0])[::-1] * np.pi / 180 

x = np.outer(np.cos(lons), np.cos(lats)).T
y = np.outer(np.sin(lons), np.cos(lats)).T
z = np.outer(np.ones(np.size(lons)), np.sin(lats)).T

ax.plot_surface(x, y, z, rstride=4, cstride=4, facecolors=bm)

# Plot Magnetic field lines. -------------------------------------------
# The magnetic field lines data needs to be converted to a rectangular
# coordinate system where the earth is centered at (0, 0, 0) and has a
# radius of 1.
#
t = np.linspace(0, 1) 
x = [0,  .5, 1,  .5,  0]
y = [0,   0, 0,   0,  0]
z = [1, 1.5, 0, -1.5, -1]

ax.plot(x, y, z, 'r-')

plt.show()
