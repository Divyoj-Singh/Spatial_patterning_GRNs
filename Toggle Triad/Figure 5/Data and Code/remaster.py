# Remastering old graphs because they had some problems

# Not to worry too much, I was a genius and had saved the data for the graphs and since I do not have to plot them, I can reuse them

# %%

# Imports
import numpy as np
import matplotlib.pyplot as plt

# First read the file using numpy
data = np.loadtxt("data.txt")
sol_p = np.array([data[:50], data[50:100], data[100:150]])

# SIMPLE PLOTS
# For each node
for i in range(3):
    # Plot the image and save it
    plt.clf()
    plt.imshow(sol_p[i], cmap="jet", interpolation="gaussian", vmin=0, vmax=1)
    plt.colorbar()
    # plt.savefig(f"Morphogen {i+1}.svg")

# Also make an RGB Plot
sol = np.moveaxis(sol_p, 0, -1)

# OVERLAY PLOT
# Make an image of out of the plot and save it
plt.clf()
plt.imshow(sol)
# plt.savefig("overlay.svg")