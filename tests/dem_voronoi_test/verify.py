import csv
import math
from pathlib import Path


def read(iteration):
    with Path(f"backup/iter={iteration}/particles.dat").open(newline="", encoding="utf-8") as stream:
        return {int(row["id"]): row for row in csv.DictReader(stream)}


initial = read(0)
final = read(100)
initial_volume = sum(float(row["total_volume"]) for row in initial.values())
final_volume = sum(float(row["total_volume"]) for row in final.values())
assert math.isclose(initial_volume, 9.0, rel_tol=1e-12, abs_tol=1e-12)
# particles.dat uses compact decimal formatting; allow its accumulated round-off.
assert math.isclose(final_volume, 9.0, rel_tol=1e-6, abs_tol=1e-5)
assert math.isclose(float(final[5]["px"]), 1.6, rel_tol=1e-10)
assert not math.isclose(float(final[5]["volume_ratio"]), 1.0, abs_tol=1e-5)
faces = Path("backup/iter=100/gas_voronoi_faces.dat")
assert faces.exists() and faces.stat().st_size > 100
print(f"PASS: conserved clipped Voronoi volume={final_volume:.9f}, moving-cell ratio={float(final[5]['volume_ratio']):.6f}")
