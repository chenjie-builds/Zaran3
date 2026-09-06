#!/usr/bin/env python3
"""Verify that a moving particle collides with an immobile DEM particle."""

import argparse
import csv
from pathlib import Path


def read_by_id(path: Path) -> dict[int, dict[str, float]]:
    result = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            key = "id" if "id" in row else "variables=id"
            result[int(row[key])] = {name: float(value) for name, value in row.items()}
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", type=Path, required=True, help="backup root")
    args = parser.parse_args()

    initial = read_by_id(args.sim / "iter=0" / "particles.dat")
    final = read_by_id(args.sim / "iter=5000" / "particles.dat")
    moving = final[0]
    fixed = final[1]

    if moving["vx"] >= -0.5:
        raise AssertionError(f"moving particle did not rebound: vx={moving['vx']}")
    for component in ("px", "py", "pz", "vx", "vy", "vz"):
        if abs(fixed[component] - initial[1][component]) > 1e-12:
            raise AssertionError(f"fixed particle changed {component}")

    print(f"moving particle final vx: {moving['vx']:.6f} m/s")
    print("PASS: fixed particle participates in contact and remains immobile")


if __name__ == "__main__":
    main()
