#!/usr/bin/env python3
"""Validate the inert single-particle DEM bounce case from Zaran backups."""

import argparse
import csv
from pathlib import Path

RADIUS = 0.01
INITIAL_Y = 0.20
RESTITUTION = 0.90
REFERENCE_PEAK = RADIUS + RESTITUTION**2 * (INITIAL_Y - RADIUS)


def read_particle(path: Path) -> tuple[float, float]:
    with path.open(newline="", encoding="utf-8") as stream:
        row = next(csv.DictReader(stream), None)
    if row is None:
        raise RuntimeError(f"empty particle backup: {path}")
    return float(row["py"]), float(row["vy"])


def load_history(backup: Path) -> list[tuple[int, float, float]]:
    history = []
    for folder in backup.glob("iter=*"):
        try:
            iteration = int(folder.name.split("=", 1)[1])
        except ValueError:
            continue
        particle_file = folder / "particles.dat"
        if particle_file.exists():
            y, vy = read_particle(particle_file)
            history.append((iteration, y, vy))
    return sorted(history)


def validate(history: list[tuple[int, float, float]]) -> None:
    if len(history) < 3:
        raise AssertionError("at least three backup states are required")

    if not any(vy < -0.1 for _, _, vy in history):
        raise AssertionError("no descending state was recorded")
    ascending_indices = [i for i, (_, _, vy) in enumerate(history) if vy > 0.1]
    if not ascending_indices:
        raise AssertionError("the particle did not rebound")

    post_bounce = history[ascending_indices[0]:]
    numerical_peak = max(y for _, y, _ in post_bounce)
    error = abs(numerical_peak - REFERENCE_PEAK)
    tolerance = 0.01

    print(f"states          : {len(history)}")
    print(f"reference peak  : {REFERENCE_PEAK:.8f} m")
    print(f"numerical peak  : {numerical_peak:.8f} m")
    print(f"absolute error  : {error:.3e} m")
    if error > tolerance:
        raise AssertionError(f"rebound peak error exceeds {tolerance} m")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", type=Path, required=True, help="backup root")
    args = parser.parse_args()
    validate(load_history(args.sim))
    print("PASS: inert DEM bounce response is physically consistent")


if __name__ == "__main__":
    main()
