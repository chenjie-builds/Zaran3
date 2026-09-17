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

# --- "每个粒子外部相切"约定 ---
# 试件半径 = 0.5 * 晶格间距 → 相邻试件粒子恰好外切（中心距 = r_a + r_b = 2r）；
# 压板球心间距恒为 2 * 压板半径 → 压板球彼此恰好外切。
radius = {int(row["id"]): float(row["radius"]) for row in particles}
position = {int(row["id"]): (float(row["px"]), float(row["py"]), float(row["pz"]))
            for row in particles}
worst = 0.0
for bond in bonds:
    a, b = int(bond["particle_a_id"]), int(bond["particle_b_id"])
    worst = max(worst, abs(math.dist(position[a], position[b]) - (radius[a] + radius[b])))
assert worst < 1e-15, f"adjacent specimen particles must be exactly tangent, worst gap {worst:.3e}"
assert abs(float(bonds[0]["rest_length"]) - 2 * float(disk[0]["radius"])) < 1e-15

platen_x = sorted(float(row["px"]) for row in top)
platen_gaps = [platen_x[i + 1] - platen_x[i] for i in range(len(platen_x) - 1)]
platen_radius = float(top[0]["radius"])
assert all(abs(gap - 2 * platen_radius) < 1e-15 for gap in platen_gaps), \
    "platen spheres must be exactly externally tangent"

# 离散球边界的封缝判据：contact_rebound_ratio 必须大于 R_b/(R_b + R_p)
specimen_radius = float(disk[0]["radius"])
criterion = platen_radius / (platen_radius + specimen_radius)
assert 0.6 > criterion, (f"the default overlap ratio 0.6 must seal the boundary "
                         f"(criterion {criterion:.4f})")

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
