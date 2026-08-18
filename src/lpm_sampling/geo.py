"""Geographic helpers, so a global frame is sampled on a sphere and not on a rectangle.

Latitude and longitude are not a metric space. Two failures matter enough to
have their own tests:

* **The antimeridian.** Longitude is circular but cut at +/-180, and Euclidean
  distance treats the cut as a real gap. Suva (+178.4) and Apia (-171.8) are
  about 1200 km apart; in degrees they look 350 apart, i.e. maximally distant.
  A sampler working in degrees therefore believes the Pacific rim is empty on
  both sides and will happily select both.
* **Anisotropy.** A degree of latitude is ~111 km everywhere, a degree of
  longitude is ~111*cos(lat) km: 111 km at the equator, 48 km at Reykjavik.
  Treating them as equivalent over-weights east-west separation with latitude.

Projecting to a plane does not fix either problem. Converting to the unit
sphere does: chord distance there is ``2*sin(great_circle/2)``, a strictly
increasing function of great-circle distance, so any algorithm that only
compares distances behaves as if it worked in kilometres, with no seam and no
pole degeneracy.
"""

from __future__ import annotations

import numpy as np

__all__ = ["EARTH_RADIUS_KM", "unit_sphere_xyz", "great_circle_km", "nearest_distance_km"]

EARTH_RADIUS_KM = 6371.0088


def unit_sphere_xyz(lat_deg, lon_deg):
    """Latitude/longitude in degrees -> ``(N, 3)`` points on the unit sphere.

    Pass the result as the coordinates for :func:`lpm_sampling.lpm.lpm`.
    """
    lat = np.radians(np.asarray(lat_deg, dtype=np.float64)).ravel()
    lon = np.radians(np.asarray(lon_deg, dtype=np.float64)).ravel()

    if lat.shape != lon.shape:
        raise ValueError("lat and lon must have the same length")

    cos_lat = np.cos(lat)

    return np.stack([cos_lat * np.cos(lon), cos_lat * np.sin(lon), np.sin(lat)], axis=1)


def great_circle_km(a_xyz, b_xyz):
    """Great-circle distance in km between unit-sphere points, elementwise over rows."""
    dot = np.clip(np.sum(np.asarray(a_xyz) * np.asarray(b_xyz), axis=-1), -1.0, 1.0)

    return EARTH_RADIUS_KM * np.arccos(dot)


def nearest_distance_km(frame_xyz, selected_xyz, chunk=2048):
    """For every frame point, great-circle km to the nearest selected point.

    Chunked so a frame of tens of thousands of points does not allocate an
    N x n matrix in one go.
    """
    frame_xyz = np.asarray(frame_xyz, dtype=np.float64)
    selected_xyz = np.asarray(selected_xyz, dtype=np.float64)
    out = np.empty(len(frame_xyz), dtype=np.float64)

    for start in range(0, len(frame_xyz), chunk):
        block = frame_xyz[start:start + chunk]
        dot = np.clip(block @ selected_xyz.T, -1.0, 1.0)
        out[start:start + len(block)] = EARTH_RADIUS_KM * np.arccos(dot.max(axis=1))

    return out
