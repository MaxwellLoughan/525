
# -*- coding: utf-8 -*-
"""
ATTENTION: While I believe that my methods for 
modeling the particle orbits are correct, I have made some 
mistake somewhere. My orbits are totally wrong, too symmetric, and 
nothing like we should expect. This document should have been uploaded
to CAVAS as a link to a running file, where I a currently debugging 
this error as of (2/25/2025).

@author: Maxwell Loughan
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from sympy import symbols, Eq, solve
import matplotlib.ticker as ticker
import matplotlib

# Constants
B0   = 7.0              # [T] reference toroidal field
J    = 20e6             # [A/m^2] current density
mu0  = 4.0*np.pi*1e-7   # [T·m/A] vacuum permeability
# Particle
m    = 1.67e-27             # [kg] proton mass
q    = 1.602e-19            # [C] proton charge
E    = 100                  # [keV] particle energy
E_J  = E*1.602e-16          # [J] particle energy in Jouls
# Particle Movement
alpha_deg = 80.                         # [deg] Pitch angle Alpha in degrees
alpha_rad = np.radians(alpha_deg)       # [rad] Pitch angle Alpha in radians
v         = np.sqrt(2 * E_J / m)      # [m/s] Total Particle Velocity
v_par     = v * np.cos(alpha_rad)       # [m/s] Parallel Velocity component
v_per     = v * np.sin(alpha_rad)       # [m/s] Perpendicular Velocity component
# Tokamak parameter [m]
R0      = 0.7        # Major radius
a       = 0.2         # Minor radius
kappa   = 1.5     # Elongation
delta   = 0.3     # Triangularity

#  Make the R, Z grid

Rv = np.linspace(0.2, 1, 400)
Zv = np.linspace(-0.3, 0.3, 400)
Rm, Zm = np.meshgrid(Rv, Zv)

# print("Rm:\n", Rm)
# print("Zm:\n", Zm)


# Compute Flux Function Psi(R,Z)
epsilon = 1e-10  # Avoid division by zero
Psi = Rm * ((mu0 * J * R0 / (4.0 * (Rm + epsilon))) * ((Rm - R0)**2 + Zm**2))

# Plot the Flux Contours
plt.figure(figsize=(7,6))
contour = plt.contourf(Rm, Zm, Psi, levels=50, cmap='viridis')
plt.colorbar(label='Flux Function (Psi)')
plt.xlabel('R (m)')
plt.ylabel('Z (m)')
plt.title('Flux Function Contours')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
############

# Now lets get our p constant for momentum

# This relies on initial conditions of a starting point (R,Z)

# ---------------------------------
# Compute Magnetic Field from ∇ × A
# ---------------------------------
B_R = (mu0 * J * R0 * Zm) / (2 * Rm**2)
B_Z = -((mu0 * J * R0) / (2 * Rm)) * (Rm - R0)
B_phi = (B0 * R0) / Rm  # Toroidal component

# Compute Total B Magnitude
B_tot = np.sqrt(B_R**2 + B_Z**2 + B_phi**2) # I will be substituting "B" in equations for B_p below 
B = np.sqrt(B_R**2 + B_Z**2) # this is technically B_p
                                  # is the component of B in the phi dir.

mu = (m*v_per**2) / (2*B)
# E1 = 0.5*m*v**2
lamb = mu/E_J

# Now we can set some inititla conditions for orbits
# remember that p_phi(R_init,Z_init) is a constant that will be for a respective orbit
# lets do 5 orbits and then plot them

# def p_phi(R_init,Z_init):
#     return R_init*m*v*(B_phi / B)*np.sqrt(1 - (B)*(lamb)) + q*Psi

# p_phi = R_init*m*v*(B_phi / B)*np.sqrt(1 - (B)*(lamb)) + q*Psi

# F = ((p_phi - q*Psi)**2)-(((m*v*R)*(B_phi/B))**2)*(1-(lamb)*B)

# F = ((p_phi - q*Rm*((mu0 * J * R0 / (4.0 * Rm)) * ((Rm - R0)**2 + Zm**2)))**2)-(((m*v*Rm)*(B_phi/B))**2)*(1-(lamb)*B)

####
initial_conditions = [ # (R_init, Z_init)
    (0.5, -0.1),  
    (0.75, 0.0),  
    (0.72, -0.05),  
    (0.78, 0.05),  
    (0.80, 0.0)   
]

plt.figure(figsize=(8,6))
color_list = ['red', 'blue', 'green', 'purple', 'orange']

for i, (R_init, Z_init) in enumerate(initial_conditions):
    # Compute p_phi at initial (R, Z)
    B_phi_init = (B0 * R0) / R_init
    Psi_init = R_init * ((mu0 * J * R0 / (4.0 * R_init)) * ((R_init - R0)**2 + Z_init**2))
    p_phi_init = R_init * m * v * (B_phi_init / B_phi) * np.sqrt(1 - B * lamb) + q * Psi_init
    # p_phi_init = R_init*m*v*(B_phi / B)*np.sqrt(1 - (B)*(lamb)) + q*Psi
    print(p_phi_init)
    # Compute F(R,Z) for this orbit
    F = ((p_phi_init - q * Psi)**2) - (((m * v * Rm) * (B_phi / B))**2 * (1 - lamb * B))

    # Plot F=0 contour
    cs = plt.contour(Rm, Zm, F, levels=[0], colors=[color_list[i % len(color_list)]], linewidths=1)
    # if len(cs.allsegs[0]) > 0:
    #     cs.collections[0].set_label(f"Orbit {i+1}: R={R_init}, Z={Z_init}")

plt.xlabel("R (m)")
plt.ylabel("Z (m)")
plt.title("F(R,Z)=0 Orbit Contours")
plt.legend()
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True)

plt.show()