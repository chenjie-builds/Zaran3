#!/usr/bin/env python3
"""Validate the 10k Brazilian disk input and optional final output."""

import csv
import math
from pathlib import Path

root = Path(__file__).resolve().parent
with (root / "particles.csv").open(newline="", encoding="utf-8") as stream:
    particles = list(csv.DictReader(stream))
with (root / "bonds.csv").open(newline="", encoding="utf-8") as stream:
    bonds = list(csv.DictReader(stream))

disk = [row for row in particles if int(row["group"]) == 0]
top = [row for row in particles if int(row["group"]) == 1]
bottom = [row for row in particles if int(row["group"]) == 2]
assert len(disk) == 10_000
assert len(top) == len(bottom) > 100
assert 29_000 < len(bonds) < 30_000
assert math.isclose(sum(float(row["px"]) for row in disk), 0.0, abs_tol=1e-12)
assert math.isclose(sum(float(row["py"]) for row in disk), 0.0, abs_tol=1e-12)
assert all(float(row["vy"]) < 0.0 for row in top)
assert all(float(row["vy"]) > 0.0 for row in bottom)

final_path = root / "backup" / "iter=10000" / "particles.dat"
if final_path.exists():
    with final_path.open(newline="", encoding="utf-8") as stream:
        final = {int(row["id"]): row for row in csv.DictReader(stream)}
    assert len(final) == len(particles)
    assert all(math.isfinite(float(row["px"])) and math.isfinite(float(row["py"]))
               for row in final.values())
    top_force = sum(float(final[int(row["id"])]["fy"]) for row in top)
    bottom_force = sum(float(final[int(row["id"])]["fy"]) for row in bottom)
    print(f"final platen reactions: top={top_force:.6e} N, bottom={bottom_force:.6e} N")
print(f"PASS: 10000 disk nodes, {len(top) + len(bottom)} platen particles, {len(bonds)} bonds")
