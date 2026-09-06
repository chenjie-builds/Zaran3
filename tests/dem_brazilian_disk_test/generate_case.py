#!/usr/bin/env python3
"""Generate an inert bonded-particle Brazilian disk and two kinematic platens."""

import argparse
import csv
import math
from pathlib import Path

DISK_RADIUS = 0.05
SPACING = 0.01
PARTICLE_RADIUS = 0.0035
PARTICLE_MASS = 0.01
PLATEN_RADIUS = 0.004
PLATEN_SPEED = 0.01


def disk_points() -> list[tuple[float, float]]:
    points = []
    row_height = math.sqrt(3.0) * 0.5 * SPACING
    rows = math.ceil(DISK_RADIUS / row_height)
    columns = math.ceil(DISK_RADIUS / SPACING) + 1
    for row in range(-rows, rows + 1):
        y = row * row_height
        offset = 0.5 * SPACING if row % 2 else 0.0
        for column in range(-columns, columns + 1):
            x = column * SPACING + offset
            if math.hypot(x, y) <= DISK_RADIUS - 0.35 * SPACING:
                points.append((x, y))
    return points


def write_particles(output: Path, points: list[tuple[float, float]]) -> None:
    top = max(y for _, y in points) + PARTICLE_RADIUS + PLATEN_RADIUS + 1e-4
    bottom = -top
    platen_count = math.ceil(DISK_RADIUS / PLATEN_RADIUS)
    platen_x = [i * 2.0 * PLATEN_RADIUS for i in range(-platen_count, platen_count + 1)]

    with (output / "particles.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["id", "group", "radius", "mass", "px", "py", "pz",
                         "vx", "vy", "vz", "ox", "oy", "oz", "motion_type"])
        for particle_id, (x, y) in enumerate(points):
            writer.writerow([particle_id, 0, PARTICLE_RADIUS, PARTICLE_MASS,
                             x, y, 0, 0, 0, 0, 0, 0, 0, 1])
        next_id = len(points)
        for x in platen_x:
            writer.writerow([next_id, 1, PLATEN_RADIUS, 1.0,
                             x, top, 0, 0, -PLATEN_SPEED, 0, 0, 0, 0, 2])
            next_id += 1
        for x in platen_x:
            writer.writerow([next_id, 2, PLATEN_RADIUS, 1.0,
                             x, bottom, 0, 0, PLATEN_SPEED, 0, 0, 0, 0, 2])
            next_id += 1


def write_bonds(output: Path, points: list[tuple[float, float]]) -> None:
    bonds = []
    for i, (xi, yi) in enumerate(points):
        for j in range(i + 1, len(points)):
            xj, yj = points[j]
            distance = math.hypot(xj - xi, yj - yi)
            if 0.98 * SPACING <= distance <= 1.02 * SPACING:
                bonds.append((len(bonds), i, j, distance, 2.0e4, 5.0e3, 2.0, 0.5, 1))
    with (output / "bonds.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["id", "particle_a_id", "particle_b_id", "rest_length",
                         "kn", "kt", "cn", "ct", "active"])
        writer.writerows(bonds)
    print(f"generated {len(points)} disk particles and {len(bonds)} bonds")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    points = disk_points()
    write_particles(args.output, points)
    write_bonds(args.output, points)


if __name__ == "__main__":
    main()
