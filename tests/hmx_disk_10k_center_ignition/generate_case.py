#!/usr/bin/env python3
"""Generate an exactly 10,000-node HMX disk on a symmetric triangular lattice."""

import csv
import math
from pathlib import Path

TARGET_NODES = 10_000
DISK_RADIUS = 0.020
THICKNESS = 0.10
DENSITY = 1890.0
YOUNG_MODULUS = 28.764e9
POISSON_RATIO = 0.20
SPECIFIC_HEAT = 1004.62
THERMAL_CONDUCTIVITY = 0.5358
REACTION_HEAT = 5.53e6
ARRHENIUS_PREFACTOR = 4.78e12
ACTIVATION_TEMPERATURE = 143.9e3 / 8.31446261815324


def make_lattice():
    # The half-cell x shift removes a point at the origin. Every lattice point
    # then has the exact inversion partner (-q-1,-r), allowing an even 10,000
    # points without breaking disk symmetry.
    limit = 90
    coordinates = [(q, r) for r in range(-limit, limit + 1)
                   for q in range(-limit, limit + 1)]
    coordinate_set = set(coordinates)
    pairs = []
    visited = set()
    for coordinate in coordinates:
        if coordinate in visited:
            continue
        q, r = coordinate
        partner = (-q - 1, -r)
        if partner not in coordinate_set:
            continue
        visited.add(coordinate)
        visited.add(partner)
        x = q + 0.5 * r + 0.5
        y = 0.5 * math.sqrt(3.0) * r
        pairs.append((x * x + y * y, coordinate, partner))
    pairs.sort(key=lambda item: (item[0], item[1]))
    selected = [coordinate for _, a, b in pairs[: TARGET_NODES // 2]
                for coordinate in (a, b)]
    max_radius_units = max(math.hypot(q + 0.5 * r + 0.5,
                                     0.5 * math.sqrt(3.0) * r)
                           for q, r in selected)
    spacing = DISK_RADIUS / (max_radius_units + 0.45)
    points = [(q, r, (q + 0.5 * r + 0.5) * spacing,
               0.5 * math.sqrt(3.0) * r * spacing) for q, r in selected]
    points.sort(key=lambda item: (item[3], item[2]))
    return points, spacing


def main():
    root = Path(__file__).resolve().parent
    points, spacing = make_lattice()
    cell_area = 0.5 * math.sqrt(3.0) * spacing * spacing
    mass = DENSITY * cell_area * THICKNESS
    interface_area = spacing * THICKNESS

    particle_header = [
        "id", "group", "radius", "mass", "px", "py", "pz",
        "vx", "vy", "vz", "ox", "oy", "oz", "motion_type",
        "young_modulus", "poisson_ratio", "friction_coeff", "restitution_coeff",
        "energetic", "temperature", "reaction_progress", "specific_heat",
        "thermal_conductivity", "reaction_heat", "arrhenius_prefactor",
        "activation_temperature", "heat_source",
    ]
    ids = {}
    with (root / "particles.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(particle_header)
        for particle_id, (q, r, x, y) in enumerate(points, 1):
            ids[(q, r)] = particle_id
            writer.writerow([
                particle_id, 1, 0.45 * spacing, mass, x, y, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1,
                YOUNG_MODULUS, POISSON_RATIO, 0.10, 0.90,
                1, 300.0, 0.0, SPECIFIC_HEAT, THERMAL_CONDUCTIVITY,
                REACTION_HEAT, ARRHENIUS_PREFACTOR, ACTIVATION_TEMPERATURE, 0.0,
            ])

    neighbor_offsets = ((1, 0), (0, 1), (-1, 1))
    bonds = []
    for q, r, _, _ in points:
        for dq, dr in neighbor_offsets:
            neighbor = (q + dq, r + dr)
            if neighbor in ids:
                bonds.append((ids[(q, r)], ids[neighbor]))
    with (root / "bonds.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["id", "particle_a_id", "particle_b_id", "rest_length",
                         "kn", "kt", "cn", "ct", "active", "conduction_area"])
        for bond_id, (a, b) in enumerate(bonds, 1):
            writer.writerow([bond_id, a, b, spacing, 1.0e7, 4.0e6,
                             20.0, 10.0, 1, interface_area])
    print(f"generated {len(points)} HMX nodes, {len(bonds)} bonds, spacing={spacing:.9e} m")


if __name__ == "__main__":
    main()
