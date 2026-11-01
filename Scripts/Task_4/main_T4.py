import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from do_time_step_T4 import do_time_step_T4

# Simulation parameters
N = 128
dx = 1.0
dt = 0.001
M = 1.0
kappa = 1.0

steps = 10000
plot_every = 2500

# Select free energy potential
potential = "triple_well"
parameter = 1.0

# Initial concentration field for Flory-Huggins
np.random.seed(42)

c = 0.5 + 0.02 * np.random.randn(N, N)

# Add small random fluctuations
c += 0.02 * np.random.randn(N, N)

# Create folder for snapshots
output_dir = Path(__file__).resolve().parent / "snapshots"
output_dir.mkdir(exist_ok=True)

# Run simulation
for step in range(steps + 1):

    if step % plot_every == 0:
        print(f"\nStep {step}")
        print(f"Minimum concentration: {np.min(c)}")
        print(f"Maximum concentration: {np.max(c)}")
        print(f"Mean concentration: {np.mean(c)}")
        print(f"Contains NaN: {np.isnan(c).any()}")
        plt.figure(figsize=(6, 5))
        plt.imshow(c, cmap="RdBu", origin="lower",
                   vmin=-1, vmax=1)
        plt.colorbar(label="Concentration")
        plt.title(f"{potential} - Step {step}")
        plt.tight_layout()

        filename = output_dir / f"{potential}_step_{step}.png"
        plt.savefig(filename, dpi=150)
        plt.close()

        print(f"Saved: {filename.name}")

    if step < steps:
        do_time_step_T4(c, dx, dt, M, kappa,potential, parameter)

        if not np.all(np.isfinite(c)):
            raise RuntimeError(f"Simulation became unstable at step {step + 1}")

print("Simulation finished!")
