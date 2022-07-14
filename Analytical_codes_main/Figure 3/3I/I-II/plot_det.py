import numpy as np
import matplotlib.pyplot as plt

x, y = np.loadtxt('TSSA_tristable_det_1D.dat', unpack=True, usecols=[0,1])

plt.plot(x, y, label="1D diffusion", color='blue', linestyle='solid', linewidth = 4, marker='o', markerfacecolor='blue', markeredgecolor="navy", markersize=12)

x, y = np.loadtxt('TSSA_tristable_det_2D.dat', unpack=True, usecols=[0,1])

plt.plot(x, y, label="2D diffusion", color='brown', linestyle='solid', linewidth = 4, marker='^', markerfacecolor='brown', markeredgecolor="maroon", markersize=12)

x, y = np.loadtxt('TSSA_tristable_det_0D.dat', unpack=True, usecols=[0,1])

plt.plot(x, y, marker="*", markersize=15, markeredgecolor="black", markerfacecolor="black", label="No diffusion")

plt.plot(x, y, marker="*", markersize=15, markeredgecolor="black", markerfacecolor="black")

plt.xlabel('k', fontsize=18, weight='bold')
plt.ylabel('Determinant', fontsize=18, weight='bold')

plt.xticks(fontsize=15, weight='bold')
plt.yticks(fontsize=15, weight='bold')

plt.rc('font', size=18)
plt.rc('font', weight='bold')
plt.legend(loc='best')
plt.tight_layout()

plt.savefig('TSSA_tri_det.png', dpi=300)
