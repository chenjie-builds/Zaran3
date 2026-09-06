#!/usr/bin/env python3
"""Check the elastic response and symmetry of the inert Brazilian disk case."""

import argparse
import csv
import math
from pathlib import Path


def particles(path: Path) -> dict[int, dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return {int(row["id"]): {key: float(value) for key, value in row.items()}
                for row in csv.DictReader(stream)}


def bonds(path: Path) -> list[tuple[int, int]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return [(int(row["particle_a_id"]), int(row["particle_b_id"]))
                for row in csv.DictReader(stream)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", type=Path, required=True)
    parser.add_argument("--case", type=Path, required=True)
    args = parser.parse_args()

    initial = particles(args.sim / "iter=0" / "particles.dat")
    final = particles(args.sim / "iter=5000" / "particles.dat")
    topology = bonds(args.case / "bonds.csv")
    disk_ids = [key for key, value in initial.items() if value["group"] == 0]
    top_ids = [key for key, value in initial.items() if value["group"] == 1]
    bottom_ids = [key for key, value in initial.items() if value["group"] == 2]

    expected_motion = 0.01 * 0.05
    top_motion = sum(initial[i]["py"] - final[i]["py"] for i in top_ids) / len(top_ids)
    bottom_motion = sum(final[i]["py"] - initial[i]["py"] for i in bottom_ids) / len(bottom_ids)
    if abs(top_motion - expected_motion) > 1e-8 or abs(bottom_motion - expected_motion) > 1e-8:
        raise AssertionError("kinematic platen displacement is incorrect")

    center_y = sum(final[i]["py"] for i in disk_ids) / len(disk_ids)
    if abs(center_y) > 2e-5:
        raise AssertionError(f"disk lost vertical symmetry: center_y={center_y}")

    initial_height = max(initial[i]["py"] for i in disk_ids) - min(initial[i]["py"] for i in disk_ids)
    final_height = max(final[i]["py"] for i in disk_ids) - min(final[i]["py"] for i in disk_ids)
    if final_height >= initial_height:
        raise AssertionError("disk did not contract along the loading diameter")

    horizontal_strains = []
    for a, b in topology:
        xa, ya = initial[a]["px"], initial[a]["py"]
        xb, yb = initial[b]["px"], initial[b]["py"]
        if math.hypot(0.5 * (xa + xb), 0.5 * (ya + yb)) > 0.015:
            continue
        if abs(yb - ya) > 0.2 * abs(xb - xa):
            continue
        l0 = math.hypot(xb - xa, yb - ya)
        l1 = math.hypot(final[b]["px"] - final[a]["px"],
                        final[b]["py"] - final[a]["py"])
        horizontal_strains.append((l1 - l0) / l0)
    mean_horizontal_strain = sum(horizontal_strains) / len(horizontal_strains)
    if mean_horizontal_strain <= 0.0:
        raise AssertionError("center region has no transverse tensile strain")

    top_reaction = sum(final[i]["fy"] for i in top_ids)
    bottom_reaction = sum(final[i]["fy"] for i in bottom_ids)
    if top_reaction <= 0.0 or bottom_reaction >= 0.0:
        raise AssertionError("platen reaction directions are incorrect")
    imbalance = abs(top_reaction + bottom_reaction) / max(abs(top_reaction), abs(bottom_reaction))
    if imbalance > 0.1:
        raise AssertionError(f"platen reactions are imbalanced: {imbalance}")

    print(f"disk particles             : {len(disk_ids)}")
    print(f"vertical diameter change   : {final_height - initial_height:.6e} m")
    print(f"center horizontal strain   : {mean_horizontal_strain:.6e}")
    print(f"top/bottom reactions       : {top_reaction:.6e}, {bottom_reaction:.6e} N")
    print(f"relative reaction imbalance: {imbalance:.3e}")
    print("PASS: elastic Brazilian disk response is symmetric and mechanically consistent")


if __name__ == "__main__":
    main()
