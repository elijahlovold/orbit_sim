import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
c = 299792458  # Speed of light (m/s)
M = 5.972e24  # Mass of the black hole (kg)
R_s = 2 * G * M / c**2  # Schwarzschild radius

def sch_radial(t, y, M):
    r, phi, pr, pphi = y
    
    # Hamiltonian equations of motion
    dr_dt = pr
    dphi_dt = pphi / r**2
    
    dpr_dt = pphi**2 / r**3 - (G * M / r**2)
    dpphi_dt = 0  # Conserved angular momentum
    
    return [dr_dt, dphi_dt, dpr_dt, dpphi_dt]

# Initial conditions: [r, phi, pr, pphi]
r0 = 20 * R_s  # Initial distance from the black hole
phi0 = 0  # Initial angle
pr0 = 0  # Initial radial momentum
pphi0 = 4.5e4  # Initial angular momentum (arbitrary)

y0 = [r0, phi0, pr0, pphi0]
t_span = (0, 1e5)
t_eval = np.linspace(*t_span, 10000)

sol = solve_ivp(sch_radial, t_span, y0, t_eval=t_eval, args=(M,), method='RK45')

# Extract solutions
r = sol.y[0]
phi = sol.y[1]

# Convert polar to Cartesian coordinates
x = r * np.cos(phi)
y = r * np.sin(phi)

# Plotting the orbit
plt.figure(figsize=(8, 8))
plt.plot(x, y, label='Star Orbit')
plt.scatter(0, 0, color='black', label='Black Hole')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.legend()
plt.title('Orbit of a Star around a Black Hole')
plt.grid(True)
plt.axis('equal')
plt.show()
