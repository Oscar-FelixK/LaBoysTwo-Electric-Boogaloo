from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from Task_1.do_time_step import do_time_step

# ----------------------------------------------------
# Parameters
# ----------------------------------------------------
N = 128

dx = 1.0
dt = 0.01

M = 1.0
kappa = 1.0

steps = 50000
sample_every = 10

# ----------------------------------------------------
# Characteristic length measurement
# ----------------------------------------------------
def defining_length(grid, spacing=1.0):
    """Estimate characteristic pattern length from the number of sign changes."""
    if grid.ndim != 2:
        raise ValueError("Expected a 2D concentration field.")

    eps = 1e-12
    sign_grid = np.signbit(grid)
    sign_grid[np.abs(grid) < eps] = False

    vertical_changes = np.count_nonzero(np.diff(sign_grid, axis=1) != 0, axis=1)
    horizontal_changes = np.count_nonzero(np.diff(sign_grid, axis=0) != 0, axis=0)

    valid_vertical = vertical_changes[vertical_changes > 0]
    valid_horizontal = horizontal_changes[horizontal_changes > 0]
    if valid_vertical.size == 0 or valid_horizontal.size == 0:
        return np.nan

    mean_vertical = valid_vertical.mean()
    mean_horizontal = valid_horizontal.mean()
    if mean_vertical == 0 or mean_horizontal == 0:
        return np.nan

    return 0.5 * spacing * (grid.shape[1] / mean_vertical + grid.shape[0] / mean_horizontal)


def fit_growth_exponent(times, lengths, min_time=100.0):
    """Fit lengths to L(t) = A * t**alpha for late times t > min_time."""
    t = np.asarray(times, dtype=float)
    L = np.asarray(lengths, dtype=float)
    mask = (t > min_time) & np.isfinite(L) & (L > 0)

    if np.count_nonzero(mask) < 2:
        return np.nan, np.nan

    slope, intercept = np.polyfit(np.log(t[mask]), np.log(L[mask]), 1)
    return slope, np.exp(intercept)


# ----------------------------------------------------
# Time integration
# ----------------------------------------------------
def run_simulation():
    c = 0.05 * np.random.randn(N, N)
    times = []
    lengths = []

    for step in range(steps + 1):
        if step > 0:
            do_time_step(c, dx, dt, M, kappa)
        if step % sample_every == 0:
            times.append(step * dt)
            lengths.append(defining_length(c, dx))

    return np.asarray(times), np.asarray(lengths)


if __name__ == "__main__":
    times, lengths = run_simulation()
    alpha, prefactor = fit_growth_exponent(times, lengths, min_time=100.0)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(times, lengths, "o-", label="measured")

    late_mask = (times > 100.0) & np.isfinite(lengths) & (lengths > 0)
    if np.count_nonzero(late_mask) >= 2 and np.isfinite(alpha) and np.isfinite(prefactor):
        fit_times = np.geomspace(times[late_mask].min(), times[late_mask].max(), 200)
        fit_lengths = prefactor * fit_times ** alpha
        ax.loglog(fit_times, fit_lengths, "--", label=fr"fit: $L(t)= {prefactor:.3g} t^{{{alpha:.3f}}}$")
    else:
        print("Warning: insufficient valid late-time data to fit L(t) = A t^alpha.")

    ax.set_xlabel("t")
    ax.set_ylabel("L(t)")
    ax.set_title("Characteristic length scale versus time")
    ax.legend()
    ax.grid(True, which="both", ls="--", alpha=0.4)

    output_path = Path(__file__).resolve().parent / "task_3_length_scale.png"
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"Saved length-scale plot to {output_path}")
    if np.isfinite(alpha) and np.isfinite(prefactor):
        print(f"Growth exponent alpha = {alpha:.6f}")
        print(f"Prefactor A = {prefactor:.6f}")

