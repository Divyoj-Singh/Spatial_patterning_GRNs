# SCSCSCSCSCSCSCSCSCSC

# CSB Lab | Indian Institute of Science, Bengaluru
# Chinmay K Haritas

# Toggle Triad ODE

# %%
print("\nSpatial GRNs - Toggle Triad in 2D \n")
print("ODE Solver")

# Import Block
import helpers
from tt_vars import *
from helpers import hill, p_i
np, pd, plt, tqdm, time, njit, prange, animator, mal, os, product = helpers.libraries()

# Parameter set details
ds_i, pset_class, databook = 3, "Triad", "data.xlsx"
r = False

# Execution timesteps
T = 200

print(f"Starting with \033[95m ({nx}, {ny}) \033[0m lattice for \033[92m {T} \033[0m units. \n")

if(r): print(f"NOTE: Coloring in repressilator mode\n")

# Setting up the system

# Computational
dx, dy, dt = 0.5, 0.5, 0.5

# Loading Parameters from Dataset
g, k, coop, thr, fold = helpers.parameter_set(ds_i, pd.read_excel(databook, sheet_name=pset_class))
print()
# Showing Dataset Properties
print(f"Dataset \033[95m{ds_i}\033[0m from \033[92m {pset_class} \033[0m")

# Modelling function
@njit()
def model(N):
    dN = np.empty(np.shape(N))
    dN[0] = (
        g[0]
        * hill(N[2], fold[2][0], coop[2][0], thr[2][0])
        * hill(N[1], fold[1][0], coop[1][0], thr[1][0])
        - k[0] * N[0]
    )
    dN[1] = (
        g[1]
        * hill(N[0], fold[0][1], coop[0][1], thr[0][1])
        * hill(N[2], fold[2][1], coop[2][1], thr[2][1])
        - k[1] * N[1]
    )
    dN[2] = (
        g[2]
        * hill(N[1], fold[1][2], coop[1][2], thr[1][2])
        * hill(N[0], fold[0][2], coop[0][2], thr[0][2])
        - k[2] * N[2]
    )
    return dN


# Initial Conditions
ics = np.random.rand(1000, 3)*(max(g/k))

# Create solution mesh
sol = np.ones((T, len(ics), nodes)) * 0.0

# Assign Initial Conditions
sol[0] = ics

# SOLVE IT!
for t in tqdm(range(1, T)):
    for j in prange(len(ics)):
        sol[t][j] = sol[t - 1][j] + dt * model(sol[t - 1][j])

# Convert to printable form
print('OKAY')
sol_p = helpers.gbyk_normalization(np.moveaxis(np.moveaxis(sol, 1, 0), 2, 0), g/k)
# sol_p = (np.moveaxis(np.moveaxis(sol, 1, 0), 2, 0))

# End states
finals_A = np.array([(sol_p[0][i][-1]) for i in range(len(ics))])
finals_B = np.array([(sol_p[1][i][-1]) for i in range(len(ics))])
finals_C = np.array([(sol_p[2][i][-1]) for i in range(len(ics))])

finals = [finals_A, finals_B, finals_C]
finalz = [finals_A-np.mean(finals_A), finals_B-np.mean(finals_B), finals_C-np.mean(finals_C)]
meanz = [np.mean(finals_A), np.mean(finals_B), np.mean(finals_C)]

# Color the plots

# Produce State Name
def state(ic):
    # Check the final state
    state = ""
    for i in range(3):
        # Added Z-Score Normalization
        if(finals[i][ic] > np.mean(finals[i])): state+="1"
        else: state+="0"
    return state

# Assign Color to State
def coloring(ic):
    st = state(ic)
    if(st == "000"): return "w"
    if(st == "001"): return "r"
    if(st == "010"): return "g"
    if(st == "011"): return "y"
    if(st == "100"): return "b"
    if(st == "101"): return "m"
    if(st == "110"): return "c"
    if(st == "111"): return "k"


# Color for each IC
colors = [coloring(i) for i in range(len(ics))]

# Plots
# Creates subplots of 2 rows and 3 columns
fig, ax = plt.subplots(2, 3, figsize=(18, 12))

# fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18,6))
for ic in tqdm(prange(len(ics))):
    # Morphogen Initial Conditions
    for i in range(nodes):
        # ax = [ax1, ax2, ax3][i]
        # Main Stuff Here - No repressilator
        
        # First row for g/k normalization
        ax[0][i].plot( ((sol_p[i][ic][0:])) , c = colors[ic])
        if not(k[i]==0):
            ax[0][i].set_title(f"Morphogen {i+1} (g/k={ (g[i]/k[i]):.2f})")
        ax[0][i].set_xlabel("T")
        ax[0][i].set_ylabel(f"Level of {i+1}")

        # Second row for z-score normalization
        ax[1][i].plot( ((sol_p[i][ic][10:]-meanz[i])) , c = colors[ic])
        ax[1][i].set_title(f"Morphogen {i+1} (Z={meanz[i]:.2f})")
        ax[1][i].set_xlabel("T")
        ax[1][i].set_ylabel(f"Level of {i+1}")

# END FORS

# Good looking plots
plt.tight_layout()
plt.show()

# We need one more plot
# Histogram of final states
fig, ax = plt.subplots(1, 1, figsize=(8, 4))
states = [state(i) for i in range(len(ics))]
ax.hist(states, bins=5)
# Title and labels
ax.set_title(f"Final State Histogram (Total = {len(ics)})")
ax.set_xlabel("State")
ax.set_ylabel("Frequency")
plt.show()

# Run thru the ICs
seen = []
for ic in tqdm(prange(len(ics))):
    # If you have seen the state before, don't do it again
    if(states[ic] in seen): continue
    # If you haven't, add it to the seen list
    else:
        seen.append(states[ic])
        print(f"a = {finals_A[ic]*(g/k)[0]:.2f}, b = {finals_B[ic]*(g/k)[1]:.2f}, c = {finals_C[ic]*(g/k)[2]:.2f}")