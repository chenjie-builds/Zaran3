import csv
import math
from pathlib import Path

with Path("backup/iter=1000/bonds.dat").open(newline="", encoding="utf-8") as stream:
    bond = next(csv.DictReader(stream))
expected = 1.0 / math.sqrt(3.0)
actual = float(bond["fracture_energy"])
assert math.isclose(actual, expected, rel_tol=1e-5), (actual, expected)
assert int(bond["active"]) == 0
assert math.isclose(float(bond["damage"]), 1.0)
assert float(bond["dissipated_fracture_energy"]) >= actual
print(f"PASS: fracture_energy={actual:.6f} J, energy criterion broke the bond")
