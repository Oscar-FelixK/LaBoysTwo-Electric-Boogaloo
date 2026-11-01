from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from calculate_mass import calc_mass
from Task_1.do_time_step import do_time_step
from calculate_free_energy import free_energy

# ----------------------------------------------------
# Parameters
# ----------------------------------------------------
N = 128              # grid size
dx = 1.0
dt = 0.01

M = 1.0              # mobility
kappa = 1.0

steps = 5000
plot_every = 500

c = 0.05*np.random.randn(N, N)

if __name__ == "__main__":
    masses = np.zeros(steps)
    free_energies = np.zeros(steps)
    for step in range(steps):
        current_mass = calc_mass(c)
        current_free_energy = np.sum(free_energy(c, kappa, dx)) * dx**2
        masses[step] = current_mass
        free_energies[step] = current_free_energy
        do_time_step(c, dx, dt, M, kappa)

    time = np.arange(steps) * dt
    fig, (mass_ax, energy_ax) = plt.subplots(
        2, 1, sharex=True, figsize=(8, 6)
    )

    mass_ax.plot(time, masses)
    mass_ax.set_title("Mass over time")
    mass_ax.set_ylabel("Mass")
    mass_ax.set_xlabel("Time")
    mass_ax.grid(True)

    energy_ax.plot(time, free_energies)
    energy_ax.set_title("Free energy over time")
    energy_ax.set_xlabel("Time")
    energy_ax.set_ylabel("Free energy")
    energy_ax.grid(True)

    fig.tight_layout()
    plt.show()