import numpy as np


def free_energy(c, kappa, dx=1.0):
    grad_x = (np.roll(c, -1, axis=0) - c) / dx
    grad_y = (np.roll(c, -1, axis=1) - c) / dx
    return (c**2 - 1) ** 2 / 4 + kappa / 2 * (grad_x**2 + grad_y**2)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from calculate_mass import calc_mass
from Task_1.laplacian import laplacian

def free_energy(c, kappa):
    free_energy = (c**2-1)**2/4 + kappa/2 * laplacian(c, 1.0)**2
    return free_energy