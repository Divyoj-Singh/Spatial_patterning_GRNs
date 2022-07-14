import numpy as np
import matplotlib.pyplot as plt

x, y = np.loadtxt('TS_regimeIII_BS1_lambda_1D.dat', unpack=True, usecols=[0,1])

plt.plot(x, y, label="1D diffusion", color='blue', linestyle='solid', linewidth = 4, marker='o', markerfacecolor='blue', markeredgecolor="navy", markersize=12)

x, y = np.loadtxt('TS_regimeIII_BS1_lambda_1D.dat', unpack=True, usecols=[0,2])

plt.plot(x, y, color='blue', linestyle='solid', linewidth = 4, marker='o', markerfacecolor='blue', markeredgecolor="navy", markersize=12)

x, y = np.loadtxt('TS_regimeIII_BS1_lambda_2D.dat', unpack=True, usecols=[0,1])

plt.plot(x, y, label="2D diffusion", color='brown', linestyle='solid', linewidth = 4, marker='^', markerfacecolor='brown', markeredgecolor="maroon", markersize=12)

x, y = np.loadtxt('TS_regimeIII_BS1_lambda_2D.dat', unpack=True, usecols=[0,2])

plt.plot(x, y, color='brown', linestyle='solid', linewidth = 4, marker='^', markerfacecolor='brown', markeredgecolor="maroon", markersize=12)

x, y = np.loadtxt('TS_regimeIII_BS1_lambda_0D.dat', unpack=True, usecols=[0,1])

plt.plot(x, y, marker="*", markersize=15, markeredgecolor="black", markerfacecolor="black", label="No diffusion")

x, y = np.loadtxt('TS_regimeIII_BS1_lambda_0D.dat', unpack=True, usecols=[0,2])

plt.plot(x, y, marker="*", markersize=15, markeredgecolor="black", markerfacecolor="black")


plt.xlabel('k', fontsize=18, weight='bold')
plt.ylabel('Eigenvalues ($\lambda_p, \lambda_m$)', fontsize=18, weight='bold')

plt.xlim([-0.5, 11])

plt.xticks(fontsize=15, weight='bold')
plt.yticks(fontsize=15, weight='bold')

plt.rc('font', size=15)
plt.rc('font', weight='bold')
plt.legend(loc='best')
plt.tight_layout()

plt.text(10.2, -0.63, '$\lambda_p$', color='blue')
plt.text(10.2, -0.73, '$\lambda_p$', color='brown')
plt.text(10.2, -1.32, '$\lambda_m$', color='blue')
plt.text(10.2, -1.42, '$\lambda_m$', color='brown')


plt.savefig('TS_regimeIII_BS1_lambda.png', dpi=300)
