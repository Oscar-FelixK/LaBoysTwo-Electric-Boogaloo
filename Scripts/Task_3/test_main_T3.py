import importlib.util
from pathlib import Path

import numpy as np

module_path = Path(__file__).resolve().parent / "main_T3.py"
spec = importlib.util.spec_from_file_location("task3_main", module_path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_defining_length_checkerboard():
    grid = np.indices((8, 8)).sum(axis=0) % 2
    grid = 2 * grid - 1

    length = module.defining_length(grid)

    assert np.isclose(length, 8 / 7.0, rtol=1e-6)


def test_fit_growth_exponent():
    times = np.array([100.0, 200.0, 400.0, 800.0])
    lengths = 2.0 * times ** 0.5

    alpha, prefactor = module.fit_growth_exponent(times, lengths, min_time=100.0)

    assert np.isclose(alpha, 0.5, rtol=1e-6)
    assert np.isclose(prefactor, 2.0, rtol=1e-6)
