'''
Updated Code as of 3/10/2025
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from scipy.interpolate import RectBivariateSpline

    # Constants
e = 1.6e-19  # positive electron charge
Ep = 100e3 * e  # 100 keV in joules
mp = 1.673e-27  # proton mass
v = np.sqrt(2 * Ep / mp)  # proton velocity
    # Magnetic constants
B0 = 7
J = 2e6
R0 = 0.7
mu0 = 4 * np.pi * 1e-7  # vacuum permeability

####### Initial conditions ########
Ri = 0.9

Zi = 0.0
###################################

# Mesh grid
Rv = np.linspace(0.3, 1.4, 200)
Zv = np.linspace(-0.5, 0.5, 201)
Rm, Zm = np.meshgrid(Rv, Zv, indexing='ij')  # indexing order

Psi = mu0 * J * R0 * ((Rm - R0)**2 + Zm**2) / (4 * Rm)

# Magnetic field components
Bp = B0 * R0 / Rm
Br = -mu0 * J * R0 * Zm / 2  # Br = -d/dZ Psi
Bz = mu0 * J * R0 * (Rm - R0) / 2  # Bz = d/dR Psi
B = np.sqrt(Bp**2 + Br**2 + Bz**2)

# Plot flux
plt.figure(1)
plt.clf()
chiv = np.arange(-0.9, 1.0, 0.2)  # range of chi values
plt.contour(Rm, Zm, Psi, 14, linestyles='-.', colors='g')
plt.xlabel('R')
plt.ylabel('Z')
plt.legend()
plt.title('Poloidal Projection of Particle Orbits')

# Interpolators for Psi, Bp, and B
interp_Psi = RectBivariateSpline(Rv, Zv, Psi)
interp_Bp = RectBivariateSpline(Rv, Zv, Bp)
interp_B = RectBivariateSpline(Rv, Zv, B)

# Loop over chi values
for chi in chiv:
    # p_phi and lambda for the initial conditions
    Psii = interp_Psi.ev(Ri, Zi)
    Bpi = interp_Bp.ev(Ri, Zi)
    Bi = interp_B.ev(Ri, Zi)
    lambdai = (1 - chi**2) / Bi  # such as (mu/E)
    Pphii = e * Psii + np.sign(chi) * mp * v * Ri * np.sqrt(1 - lambdai * Bi) * Bpi / Bi

    # orbit equation
    F = (Pphii - e * Psi)**2 - (mp * v * Rm * Bp / B)**2 * (1 - lambdai * B)
    
    # Contour plot
    if chi < 0:
        plt.contour(Rm, Zm, F, levels=[0, 1e5], colors='r')
    else:
        plt.contour(Rm, Zm, F, levels=[0, 1e5], colors='b')
        
# Legend
legend_colors = ['r', 'b', 'g']
legend_labels = ['Positive Orbits', 'Negitive Orbits', 'Flux Surfaces' ]
legend_lines = [mlines.Line2D([], [], color=color, label=label) 
                for color, label in zip(legend_colors, legend_labels)]
plt.legend(handles=legend_lines, loc="upper right")

plt.show()


