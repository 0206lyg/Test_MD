"""Seeded placement of disjoint spherical W1 pockets in a spherical wet droplet."""

import numpy as np
from scipy.optimize import minimize


def _ball_points(rng, count, radius):
    """Independent uniform-in-volume points, measured from the droplet centre."""
    u, v, w = rng.random((3, count))
    z = 2.0 * u - 1.0
    azimuth = 2.0 * np.pi * v
    radial_xy = np.sqrt(1.0 - z * z)
    directions = np.column_stack((radial_xy * np.cos(azimuth),
                                  radial_xy * np.sin(azimuth), z))
    return directions * (radius * np.cbrt(w))[:, None]


def _first_clear(candidates, others, min_distance_sq):
    if not len(others):
        return candidates[0]
    delta = candidates[:, None, :] - others[None, :, :]
    d2 = np.einsum("ijk,ijk->ij", delta, delta)
    allowed = np.flatnonzero(np.all(d2 >= min_distance_sq, axis=1))
    return candidates[allowed[0]] if len(allowed) else None


def _untangle(rng, n, centre_limit, separation):
    """Find a feasible starting configuration when sequential insertion jams."""
    left, right = np.triu_indices(n, 1)
    target = separation / centre_limit

    def penalty(flat_centres):
        xyz = flat_centres.reshape(n, 3)
        radii = np.linalg.norm(xyz, axis=1)
        delta = xyz[left] - xyz[right]
        distance = np.linalg.norm(delta, axis=1)
        # A tiny buffer is used only by this numerical initializer. Its output
        # must meet the original constraints before sampling is allowed to start.
        excess = np.maximum(radii - (1.0 - 1.e-6), 0.0)
        overlap = np.maximum(target + 1.e-6 - distance, 0.0)
        derivative = delta * (overlap / np.maximum(distance, 1.e-15))[:, None]
        grad = xyz * (excess / np.maximum(radii, 1.e-15))[:, None]
        np.add.at(grad, left, -derivative)
        np.add.at(grad, right, derivative)
        return 0.5 * (excess @ excess + overlap @ overlap), grad.ravel()

    for _ in range(12):
        guess = _ball_points(rng, n, 1.0)
        result = minimize(penalty, guess.ravel(), jac=True, method="L-BFGS-B",
                          options={"ftol": 1.e-16, "gtol": 1.e-9, "maxiter": 300})
        xyz = result.x.reshape(n, 3) * centre_limit
        if (np.all(np.linalg.norm(xyz, axis=1) <= centre_limit)
                and np.all(np.linalg.norm(xyz[left] - xyz[right], axis=1) >= separation)):
            return xyz
    return None


def place_pockets(rng, n, wet_radius_nm, pocket_radius_nm, gap_nm):
    """Return n centres in nm, relative to the wet droplet centre.

    Both the outer oil film and the gap between two W1 surfaces must be at
    least gap_nm. Uniform-volume trial positions carry no preferred radius
    or direction. Exclusion at finite volume fraction can affect the final
    radial distribution. This is an initial packing, not an emulsion model.
    """
    if not isinstance(n, (int, np.integer)) or not 1 <= n <= 10:
        raise ValueError("W1 pocket count must be an integer from 1 to 10")
    if (not np.isfinite([wet_radius_nm, pocket_radius_nm, gap_nm]).all()
            or min(wet_radius_nm, pocket_radius_nm) <= 0 or gap_nm < 0):
        raise ValueError("W1 placement needs finite positive radii and a nonnegative oil gap")

    centre_limit = wet_radius_nm - pocket_radius_nm - gap_nm
    separation = 2.0 * pocket_radius_nm + gap_nm
    if centre_limit < 0 or (n > 1 and 2.0 * centre_limit <= separation):
        raise ValueError(
            f"No random packing fits {n} W1 pockets of radius {pocket_radius_nm:.3f} nm "
            f"in a wet radius of {wet_radius_nm:.3f} nm with a {gap_nm:.3f} nm oil gap"
        )
    if n == 1:
        return _ball_points(rng, 1, centre_limit)

    # Each restart discards the whole partial packing. Adding a pocket never
    # changes its radius, the specified clearances, or the requested count.
    min_distance_sq = separation ** 2
    centres = np.empty((n, 3))
    for _ in range(16):
        for i in range(n):
            trial = _ball_points(rng, 256, centre_limit)
            point = _first_clear(trial, centres[:i], min_distance_sq)
            if point is None:
                break
            centres[i] = point
        else:
            break
    else:
        centres = _untangle(rng, n, centre_limit, separation)
        if centres is None:
            raise ValueError(
                f"Could not place {n} W1 pockets (radius {pocket_radius_nm:.3f} nm, "
                f"wet radius {wet_radius_nm:.3f} nm, oil gap {gap_nm:.3f} nm) "
                "after random restarts. Change the diameter, count, or seed."
            )

    # Mix the initial packing. Global relocations sample uniformly in the
    # available volume conditional on the other centres. When most of that
    # volume is blocked, a symmetric local proposal can still move a pocket.
    labels = np.arange(n)
    for _ in range(64):
        for i in rng.permutation(n):
            trial = _ball_points(rng, 32, centre_limit)
            others = centres[labels != i]
            point = _first_clear(trial, others, min_distance_sq)
            if point is None:
                local = centres[i] + _ball_points(rng, 1, 0.12 * centre_limit)
                if np.linalg.norm(local[0]) <= centre_limit:
                    point = _first_clear(local, others, min_distance_sq)
            if point is not None:
                centres[i] = point
    return centres
