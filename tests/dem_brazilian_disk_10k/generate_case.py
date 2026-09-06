#!/usr/bin/env python3
"""Generate an exactly 10,000-node bonded Brazilian disk and two platens."""

import csv
import math
from pathlib import Path

TARGET_NODES = 10_000
DISK_RADIUS = 0.05
THICKNESS = 0.01
DENSITY = 2500.0
PLATEN_SPEED = 0.05


def make_lattice():
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
        visited.update((coordinate, partner))
        x = q + 0.5 * r + 0.5
        y = 0.5 * math.sqrt(3.0) * r
        pairs.append((x * x + y * y, coordinate, partner))
    pairs.sort(key=lambda item: (item[0], item[1]))
    selected = [coordinate for _, a, b in pairs[: TARGET_NODES // 2]
                for coordinate in (a, b)]
    max_radius_units = max(math.hypot(q + 0.5 * r + 0.5,
                                     0.5 * math.sqrt(3.0) * r)
                           for q, r in selected)
    spacing = DISK_RADIUS / (max_radius_units + 0.40)
    points = [(q, r, (q + 0.5 * r + 0.5) * spacing,
               0.5 * math.sqrt(3.0) * r * spacing) for q, r in selected]
    points.sort(key=lambda item: (item[3], item[2]))
    return points, spacing


def main():
    root = Path(__file__).resolve().parent
    points, spacing = make_lattice()
    particle_radius = 0.35 * spacing
    platen_radius = 0.40 * spacing
    cell_area = 0.5 * math.sqrt(3.0) * spacing * spacing
    particle_mass = DENSITY * cell_area * THICKNESS
    y_edge = max(y for _, _, _, y in points)
    platen_y = y_edge + particle_radius + platen_radius + 0.02 * spacing
    platen_count = math.ceil(DISK_RADIUS / (2.0 * platen_radius))
    platen_x = [2.0 * platen_radius * i for i in range(-platen_count, platen_count + 1)
                if abs(2.0 * platen_radius * i) <= DISK_RADIUS]

    ids = {}
    with (root / "particles.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["id", "group", "radius", "mass", "px", "py", "pz",
                         "vx", "vy", "vz", "ox", "oy", "oz", "motion_type"])
        for particle_id, (q, r, x, y) in enumerate(points, 1):
            ids[(q, r)] = particle_id
            writer.writerow([particle_id, 0, particle_radius, particle_mass,
                             x, y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1])
        next_id = TARGET_NODES + 1
        for x in platen_x:
            writer.writerow([next_id, 1, platen_radius, 1.0, x, platen_y, 0.0,
                             0.0, -PLATEN_SPEED, 0.0, 0.0, 0.0, 0.0, 2])
            next_id += 1
        for x in platen_x:
            writer.writerow([next_id, 2, platen_radius, 1.0, x, -platen_y, 0.0,
                             0.0, PLATEN_SPEED, 0.0, 0.0, 0.0, 0.0, 2])
            next_id += 1

    bonds = []
    for q, r, _, _ in points:
        for dq, dr in ((1, 0), (0, 1), (-1, 1)):
            neighbor = (q + dq, r + dr)
            if neighbor in ids:
                bonds.append((ids[(q, r)], ids[neighbor]))
    with (root / "bonds.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["id", "particle_a_id", "particle_b_id", "rest_length",
                         "kn", "kt", "cn", "ct", "active"])
        for bond_id, (a, b) in enumerate(bonds, 1):
            writer.writerow([bond_id, a, b, spacing, 2.0e4, 5.0e3, 0.1, 0.025, 1])
    print(f"generated {TARGET_NODES} disk nodes, {len(platen_x) * 2} platen particles, "
          f"{len(bonds)} bonds, spacing={spacing:.9e} m")


if __name__ == "__main__":
    main()
