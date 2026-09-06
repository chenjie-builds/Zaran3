import csv
import math
from pathlib import Path


def load_particles(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


initial = load_particles(Path("backup/iter=0/particles.dat"))
final = load_particles(Path("backup/iter=8000/particles.dat"))
assert len(initial) == len(final) >= 80

def radius(row):
    return math.hypot(float(row["px"]), float(row["py"]))

center = [r for r in final if radius(r) < 0.003]
outer = [r for r in final if radius(r) > 0.012]
center_alpha = sum(float(r["reaction_progress"]) for r in center) / len(center)
outer_reacted = sum(float(r["reaction_progress"]) > 0.001 for r in outer)
max_pressure = max(float(r["gas_pressure"]) for r in final)
reacted = sum(float(r["reaction_progress"]) > 0.001 for r in final)

snapshots = sorted(Path("backup").glob("iter=*/particles.dat"))
max_core_increment = 0.0
max_neighbor_increment = 0.0
max_internal_transfer = 0.0
for snapshot in snapshots:
    for row in load_particles(snapshot):
        max_core_increment = max(max_core_increment, float(row["core_burn_increment"]))
        max_neighbor_increment = max(
            max_neighbor_increment, float(row["neighbor_burn_increment"])
        )
        max_internal_transfer = max(
            max_internal_transfer, abs(float(row["internal_heat_transfer"]))
        )

with Path("backup/iter=8000/bonds.dat").open(newline="", encoding="utf-8") as stream:
    bond_rows = list(csv.DictReader(stream))
broken = sum(int(r["active"]) == 0 for r in bond_rows)

assert center_alpha > 0.90, center_alpha
assert reacted > len(final) // 3, reacted
assert outer_reacted > 0, outer_reacted
assert max_pressure > 9.0e4, max_pressure
assert broken > 0, broken
assert max_core_increment > 0.0, max_core_increment
assert max_neighbor_increment > 0.0, max_neighbor_increment
assert max_internal_transfer > 0.0, max_internal_transfer
print(
    f"PASS: nodes={len(final)}, reacted={reacted}, center_alpha={center_alpha:.3f}, "
    f"outer_reacted={outer_reacted}, max_pressure={max_pressure:.3e} Pa, broken_bonds={broken}, "
    f"max_core_dalpha={max_core_increment:.3e}, "
    f"max_neighbor_dalpha={max_neighbor_increment:.3e}"
)
