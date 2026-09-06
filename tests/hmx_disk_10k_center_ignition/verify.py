#!/usr/bin/env python3
"""Validate the generated 10k HMX input and, when present, its final output."""

import csv
import math
from pathlib import Path

root = Path(__file__).resolve().parent
with (root / "particles.csv").open(newline="", encoding="utf-8") as stream:
    particles = list(csv.DictReader(stream))
with (root / "bonds.csv").open(newline="", encoding="utf-8") as stream:
    bonds = list(csv.DictReader(stream))

assert len(particles) == 10_000
assert 29_000 < len(bonds) < 30_000
assert all(math.isclose(float(row["reaction_heat"]), 5.53e6) for row in particles)
assert all(math.isclose(float(row["arrhenius_prefactor"]), 4.78e12) for row in particles)
assert all(math.isclose(float(row["activation_temperature"]), 17307.1918906483,
                        rel_tol=1e-12) for row in particles)
assert min(math.hypot(float(row["px"]), float(row["py"])) for row in particles) < 2.0e-4

backup_root = root / "backup_lsm_fixed"
iterations = sorted(
    ((path / "particles.dat").stat().st_mtime,
     int(path.name.split("=", 1)[1]), path)
    for path in backup_root.glob("iter=*")
    if path.is_dir() and path.name.split("=", 1)[1].isdigit()
    and (path / "particles.dat").exists()
) if backup_root.exists() else []
if iterations:
    _, iteration, final_dir = iterations[-1]
    final_path = final_dir / "particles.dat"
    with final_path.open(newline="", encoding="utf-8") as stream:
        final = list(csv.DictReader(stream))
    assert len(final) == 10_000
    numeric_fields = ("px", "py", "vx", "vy", "temperature", "reaction_progress",
                      "volume_ratio", "gas_pressure")
    assert all(math.isfinite(float(row[field])) for row in final for field in numeric_fields)
    faces = final_dir / "gas_voronoi_faces.dat"
    assert faces.exists() and faces.stat().st_size > 100
    reacted = sum(float(row["reaction_progress"]) > 0.001 for row in final)
    print(f"iter={iteration} output: reacted nodes={reacted}/10000, dynamic Voronoi faces present")
print(f"PASS: 10000 HMX nodes, {len(bonds)} bonds")
