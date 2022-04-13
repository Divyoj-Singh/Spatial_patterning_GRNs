"""# Helper Functions
Provides
1. Parameter Set Initialisation
2. Functional Definition
3. Plotting tools
"""

# SCSCSCSCSC

# CSB Lab | Indian Institute of Science, Bengaluru
# Chinmay K Haritas

# Toggle Triad in 2D with Diffusion
# Helpers Functions

# Import Block
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from numba import njit, prange
import matplotlib.animation as animator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from IPython.display import display
import seaborn as sns
from tt_vars import *
import os
from itertools import product

sns.set_theme(context="notebook", style="ticks", font="arial", font_scale=1.2)
nowString = time.strftime("_%Y_%b_%d__%H_%M_%S")

# Returns all libraries
def libraries(): 
    """ # Libraries
    Libraries for Spatial GRNs
    """
    return np, pd, plt, tqdm, time, njit, prange, animator, make_axes_locatable, os, product



# Return Paramaters
def parameter_set(i: int, ps: pd.DataFrame):
    """
    # Parameters
    Returns the interactions variables for the the given set id, sheet_name and file"""
    global g
    global k
    global coop
    global thr
    global fold

    ds_i = i
    g = np.array([ps["Prod_of_A"][ds_i], ps["Prod_of_B"][ds_i], ps["Prod_of_C"][ds_i]])
    # Degradation
    k = np.array([ps["Deg_of_A"][ds_i], ps["Deg_of_B"][ds_i], ps["Deg_of_C"][ds_i]])

    # Cooperativity
    coop = np.array(
        [
            [ps["Num_of_AToA"][ds_i], ps["Num_of_AToB"][ds_i], ps["Num_of_AToC"][ds_i]],
            [ps["Num_of_BToA"][ds_i], ps["Num_of_BToB"][ds_i], ps["Num_of_BToC"][ds_i]],
            [ps["Num_of_CToA"][ds_i], ps["Num_of_CToB"][ds_i], ps["Num_of_CToC"][ds_i]],
        ]
    )

    # Threshold
    thr = np.array(
        [
            [ps["Trd_of_AToA"][ds_i], ps["Trd_of_AToB"][ds_i], ps["Trd_of_AToC"][ds_i]],
            [ps["Trd_of_BToA"][ds_i], ps["Trd_of_BToB"][ds_i], ps["Trd_of_BToC"][ds_i]],
            [ps["Trd_of_CToA"][ds_i], ps["Trd_of_CToB"][ds_i], ps["Trd_of_CToC"][ds_i]],
        ]
    )

    # Fold
    fold = np.array(
        [
            [ps["Inh_of_AToA"][ds_i], ps["Inh_of_AToB"][ds_i], ps["Inh_of_AToC"][ds_i]],
            [ps["Inh_of_BToA"][ds_i], ps["Inh_of_BToB"][ds_i], ps["Inh_of_BToC"][ds_i]],
            [ps["Inh_of_CToA"][ds_i], ps["Inh_of_CToB"][ds_i], ps["Inh_of_CToC"][ds_i]],
        ]
    )

    # Plot the image of the parameter sets

    # print("Starter of dataset", g[0])
    return g, k, coop, thr, fold


# Boundary Conditions
@njit
def p_i(x):
    """Input a [x,y] and recieve the expected coordinates according to periodic boundary conditions"""
    # X is any point to be accessed in space

    # Neumann
    # Send to first if before first asked
    # Send to last if after last asked
    y = np.copy(x)
    if x[0] == -1: y[0] = 0  # X
    if x[0] == nx + 1: y[0] = nx  # X
    if x[1] == -1: y[1] = 0  # Y
    if x[1] == ny + 1: y[1] = ny  # Y
    return y


# Hill Function
@njit
def hill(N, fold, n, thr):
    """# Shifted Hill Function
    1. Simulates the connection between two nodes.
    2. Below 1, lambda represses and above, it activates.
    """
    n_hill = 1 / (1 + (N / thr) ** n)
    return n_hill + fold * (1 - n_hill)


# IC Patching
def patch(x1, x2, y1, y2, fill_value, r_max):
    """Enter interval of patch and a function(x,y) to compute initial conditions at (x,y). """
    sol_1 = np.zeros((nx, ny, nodes))
    for x in range(x1, x2):
        for y in range(y1, y2):
            sol_1[x][y] = fill_value(x, y, r_max)
    sol_1_p = np.moveaxis(np.moveaxis(sol_1, 0, -1), 0, -1)
    return sol_1


# Special Patchers

# Uniform Fill
def fill_unif(u):
    """Fill space with uniform value 'u' """
    return lambda x, y: np.array([u[0], u[1], u[2]])


# Random Fill
def fill_rand(x, y, r_max):
    """Fill space with random integers in the interval"""
    return np.random.rand(nodes) * 30


# Custom Fill
def fill_custom(x, y):
    """Fill space with custom function"""
    return [0, 9, 0]

# G By k Normalization
def gbyk_normalization(sol_p, gbyk):
    """Performs g/k normalization for better data representation."""
    for i in range(len(sol_p)):
        sol_p[i] /= gbyk[i]
    return sol_p

# Plotting tools
def animation(sol_p, save=False, fol=""):
    # Artists
    scrubby = 10

    fig, axis = plt.subplots(1, 1, figsize=(5, 5))
    img = axis.imshow(sol_p[0][0], cmap="jet", interpolation="gaussian")

    def lilly(s):
        probe = s * scrubby
        img.set_data((sol_p[0][probe]))
        # Add colorbar
        img.autoscale()
        return img

    anim = animator.FuncAnimation(
        fig, lilly, range(len(sol_p[0]) // scrubby), interval=10
    )
    plt.show()
    if(save):
        print(f"Saving to {fol}/RGB_{nowString}.mp4")
        anim.save(f"{fol}/animation_{nowString}.mp4", writer=animator.FFMpegWriter())
    # ANIMATION ENDS HERE

def rgb_plot(sol_p, save=False, fol = ""):
    """Plots the solution in RGB space"""
    # First convert sol_p to sol
    sol = np.moveaxis(sol_p, 0, -1)
    scrubby = 10

    fig, axis = plt.subplots(1, 1, figsize=(5, 5))
    img = axis.imshow(sol[0], cmap="jet", interpolation="gaussian")

    # Artists
    def lilly(s):
        probe = s * scrubby
        img.set_data(sol[probe])
        # img.set_title(f"At frame {probe}")
        return img

    anim = animator.FuncAnimation(
        fig, lilly, range(len(sol) // scrubby), interval=10
    )
    plt.show()
    if(save):
        print(f"Saving to {fol}/RGB_{nowString}.mp4")
        anim.save(f"{fol}/animationRGB_{nowString}.mp4", writer=animator.FFMpegWriter())
    # RGB PLOT ENDS HERE

# Probe a single point
def probe_point(sol_p, x, y, gbyk):
    """Probes a single point through time."""
    # First convert sol_p to sol

    # nodes, T, nx, ny
    # nx, ny, T, nodes
    sol = np.moveaxis(np.moveaxis(sol_p, 0, -1), 0, -1)
    
    probe_data = sol[x][y]
    print(np.shape(probe_data)) # Will be 3xT

    # Plot the three series in probe_data
    fig, ax = plt.subplots()
    ax.plot(probe_data[0]/gbyk[0], label="A")
    ax.plot(probe_data[1]/gbyk[1], label="B")
    ax.plot(probe_data[2]/gbyk[2], label="C")
    ax.legend()
    
    # Set appropriate title
    plt.title(f"At point ({x}, {y})")

    plt.show()
    # END OF PROBE_POINT

# Save a single image
def save_image(img, savepath):
    """Saves an image to a file"""
    # Wipeout figure anyways
    plt.clf()
    # Render image
    plt.imshow(img, cmap="jet", interpolation="gaussian", vmin=0, vmax=1)
    plt.colorbar()
    plt.savefig(savepath)
    # SAVE IMG ENDS HERE

# Means of life
def means(sol_p):
    return np.average(sol_p[0][-1]), np.average(sol_p[1][-1]), np.average(sol_p[2][-1])

# State identifier
def state(x, means):
    # Check the final state
    state = ""
    for i in range(3):
        # Added Z-Score Normalization
        if(x[i]-means[i] > 0):
            state += "1"
        else:
            state += "0"
    return state


# Saving Tools
def save_adv(savePath, plotName, sol_p, pset_deets, snapshots, save_anim):
    # Time of saving
    folderName = savePath + "/" + plotName + nowString

    # Create new folder with timestamp and plot name
    os.mkdir(folderName)

    # Parallely add the data
    if(snapshots):
        # Repressilator

        # Restructured solution for plotting
        sol = np.moveaxis(sol_p, 0, -1)
        # First save the initial conditions for the repressilator
        rep_i = open(folderName + "/repressilator_init.txt", "a")
        # Saving node by node
        np.savetxt(rep_i, sol_p[0][0])
        np.savetxt(rep_i, sol_p[1][0])
        np.savetxt(rep_i, sol_p[2][0])
        rep_i.close()
        # Shutter speed
        # Save every 20th image from the midpoint of the simulation
        for t in np.arange(len(sol_p[0])-2000, len(sol_p[0]), 20):
            # For each morphogen
            for m in range(3):
                save_with_name = f"{folderName}/Morphogen {m+1} @ {t}."
                save_image(sol_p[m][t], save_with_name+"svg")
                save_image(sol_p[m][t], save_with_name+"png")
            # And a overlay as well
            save_image(sol[t], f"{folderName}/Overlay @ {t}.svg")
            save_image(sol[t], f"{folderName}/Overlay @ {t}.png")
    else:
        # Saving it
        data = open((folderName+"/data.txt"), "a")
        for i in range(3):
            save_with_name = f"{folderName}/Morphogen {i+1}."
            save_image(sol_p[i][-1], save_with_name+"svg")
            save_image(sol_p[i][-1], save_with_name+"png")
            np.savetxt(data, sol_p[i][-1])
        data.close()    

    # Save the animation !!!

    # RGB Plot
    if(save_anim):
        rgb_plot(sol_p, save=True, fol = folderName)
        animation(sol_p, save=True, fol = folderName)
    
    # Synthesize the metadata file and save it
    # Metadata format, look at metadata_template.txt
    md_temp = f"""SCSCSCSCSCSCSCSCSCSC\n\nAuthor:\nChinmay K Haritas,\n@CSB Lab | Indian Institute of Science, Bengaluru\n\nDate:\n{time.strftime("%Y %b %d - %H:%M:%S")}\n\nPlot Name:\n{plotName}\n\nParameter Set:\nSheet Name: {pset_deets['sheet']}\nIndex: {pset_deets['pset_id']} (Row + 2)\nDiffusion: {pset_deets['diff_coeff']} units\n"""
    with open((folderName+"/Meta.txt"), "w") as meta:
        meta.write(md_temp)
    # SAVE ADV ENDS HERE