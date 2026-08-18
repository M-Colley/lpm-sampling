"""Spatially balanced probability sampling.

The Local Pivotal Method draws a sample that is spread over space while
honouring inclusion probabilities exactly, so the result is a genuine
probability sample with usable design weights (``1 / pi``).

    >>> import numpy as np
    >>> from lpm_sampling import lpm2, pi_from_size, unit_sphere_xyz
    >>> rng = np.random.default_rng(0)
    >>> lat, lon = rng.uniform(-60, 60, 500), rng.uniform(-180, 180, 500)
    >>> pi, certainty = pi_from_size(rng.gamma(2.0, 5e4, 500), n=25)
    >>> sample = lpm2(pi, unit_sphere_xyz(lat, lon), rng=rng)
    >>> len(sample)
    25

Reference: Grafstrom, A., Lundstrom, N.L.P. & Schelin, L. (2012), "Spatially
balanced sampling through the pivotal method", Biometrics 68(2), 514-520.
"""

from .balance import spatial_balance, systematic_pps
from .geo import EARTH_RADIUS_KM, great_circle_km, nearest_distance_km, unit_sphere_xyz
from .lpm import lpm, lpm1, lpm2
from .probabilities import pi_from_size

__version__ = "0.1.0"

__all__ = [
    "lpm",
    "lpm1",
    "lpm2",
    "pi_from_size",
    "spatial_balance",
    "systematic_pps",
    "unit_sphere_xyz",
    "great_circle_km",
    "nearest_distance_km",
    "EARTH_RADIUS_KM",
    "__version__",
]
