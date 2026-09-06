import csv
import math
from pathlib import Path


def read_particles(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return {int(row["id"]): row for row in csv.DictReader(stream)}


initial = read_particles(Path("backup/iter=0/particles.dat"))
final = read_particles(Path("backup/iter=10/particles.dat"))
initial_reactive = initial[1]
final_reactive = final[1]
alpha = float(final_reactive["reaction_progress"])
expected = 1.0 - (1.0 - 10.0e-3) ** 10
assert math.isclose(alpha, expected, rel_tol=2e-5), (alpha, expected)
assert float(final_reactive["gas_internal_energy"]) > 0.0
assert float(final_reactive["gas_temperature"]) > float(final_reactive["temperature"])
assert int(final_reactive["phase"]) == 1

hot_drop = float(initial[2]["temperature"]) - float(final[2]["temperature"])
cold_rise = float(final[3]["temperature"]) - float(initial[3]["temperature"])
assert hot_drop > 0.0 and cold_rise > 0.0
assert math.isclose(hot_drop, cold_rise, rel_tol=1e-9, abs_tol=1e-9)
print(
    f"PASS: alpha={alpha:.7f}, gas_temperature={float(final_reactive['gas_temperature']):.3f} K, "
    f"conduction_delta={cold_rise:.6f} K"
)
