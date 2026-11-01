
from pathlib import Path
import sys

# Make the Scripts directory available for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Task_1.laplacian import laplacian
from free_energy_T4 import df_dc


def do_time_step_T4(c, dx, dt, M, kappa, potential, parameter=1.0):
    
#Perform one Cahn-Hilliard time step with a selected free energy.

    # Laplacian of concentration
    lap_c = laplacian(c, dx)

    # Chemical potential using the selected free energy
    mu = df_dc(c, potential, parameter) - kappa * lap_c

    # Update concentration
    c += dt * M * laplacian(mu, dx)

    return c

#test mass conservation

"""if __name__ == "__main__":
    import numpy as np

    c = 0.05 * np.random.randn(128, 128)

    print("Initial mean concentration:", np.mean(c))

    do_time_step_T4(c, 1.0, 0.01, 1.0, 1.0, "triple_well")

    print("Mean concentration after one step:", np.mean(c))
"""