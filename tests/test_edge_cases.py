"""Boundary behaviour the reference implementation must also satisfy."""

import os
import sys

import numpy as np

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src")

if SRC not in sys.path:
    sys.path.insert(0, SRC)

from lpm_sampling import lpm, lpm1, lpm2, pi_from_size  # noqa: E402


def coords(n, seed=0):
    return np.random.default_rng(seed).uniform(0.0, 100.0, size=(n, 2))


def test_single_unit_sample():
    pi = np.full(50, 1.0 / 50.0)
    picked = np.zeros(50)

    for seed in range(2000):
        selected = lpm2(pi, coords(50), rng=np.random.default_rng(seed))
        assert len(selected) == 1
        picked[selected] += 1

    frequency = picked / 2000.0
    se = np.sqrt(pi * (1 - pi) / 2000.0)
    assert np.max(np.abs(frequency - pi) / se) < 4.5


def test_census_selects_everything():
    pi = np.ones(40)
    selected = lpm2(pi, coords(40), rng=np.random.default_rng(1))
    assert len(selected) == 40


def test_forced_in_and_forced_out_are_respected():
    pi = np.concatenate([np.ones(5), np.zeros(10), np.full(35, 10.0 / 35.0)])

    for seed in range(200):
        selected = set(lpm2(pi, coords(50), rng=np.random.default_rng(seed)).tolist())
        assert len(selected) == 15
        assert set(range(5)).issubset(selected), "a certainty unit was dropped"
        assert not selected & set(range(5, 15)), "a zero-probability unit was selected"


def test_duplicate_coordinates_do_not_hang():
    pi = np.full(60, 0.25)
    stacked = np.zeros((60, 2))
    selected = lpm2(pi, stacked, rng=np.random.default_rng(3))
    assert len(selected) == 15


def test_both_variants_agree_on_the_design():
    rng = np.random.default_rng(5)
    pi, _ = pi_from_size(rng.gamma(2.0, 1000.0, 120), 20)
    points = coords(120, seed=5)

    for method in ("lpm1", "lpm2"):
        for seed in range(50):
            assert len(lpm(pi, points, method=method,
                           rng=np.random.default_rng(seed))) == 20


def test_rejects_impossible_inputs():
    for bad_pi in (np.array([0.5, 0.7]), np.array([-0.1, 1.1])):
        try:
            lpm2(bad_pi, coords(2))
        except ValueError:
            continue
        raise AssertionError("accepted invalid pi %s" % bad_pi)


def test_seeded_runs_are_reproducible():
    pi = np.full(80, 0.25)
    points = coords(80, seed=8)
    a = lpm2(pi, points, rng=np.random.default_rng(42))
    b = lpm2(pi, points, rng=np.random.default_rng(42))
    assert np.array_equal(a, b), "same seed must give the same sample"
    assert not np.array_equal(a, lpm1(pi, points, rng=np.random.default_rng(42)))


if __name__ == "__main__":
    failures = 0

    for name, test in sorted(globals().items()):
        if name.startswith("test_") and callable(test):
            try:
                test()
                print("PASS  %s" % name)
            except AssertionError as error:
                failures += 1
                print("FAIL  %s\n      %s" % (name, error))

    raise SystemExit(1 if failures else 0)
