"""Generate a small 2-D hexagonal lattice disk for the center-ignition regression."""
import csv
import math
from pathlib import Path


SPACING = 0.004
DISK_RADIUS = 0.020
THICKNESS = 0.10
DENSITY = 1890.0
YOUNG_MODULUS = 28.764e9
POISSON_RATIO = 0.20
SPECIFIC_HEAT = 1004.62
THERMAL_CONDUCTIVITY = 0.5358
REACTION_HEAT = 5.53e6
ARRHENIUS_PREFACTOR = 4.78e12
# 输入激活能 Ea=143.9 kJ/mol，程序保存 Ta=Ea/R（K）。
ACTIVATION_TEMPERATURE = 143.9e3 / 8.31446261815324
CELL_AREA = math.sqrt(3.0) * SPACING * SPACING / 2.0
MASS = DENSITY * CELL_AREA * THICKNESS
INTERFACE_AREA = SPACING * THICKNESS


points = []
limit = math.ceil(DISK_RADIUS / SPACING) + 1
for row in range(-limit, limit + 1):
    y = row * SPACING * math.sqrt(3.0) / 2.0
    for col in range(-limit, limit + 1):
        x = (col + 0.5 * (row & 1)) * SPACING
        if x * x + y * y <= DISK_RADIUS * DISK_RADIUS + 1e-15:
            points.append((x, y))
points.sort(key=lambda p: (p[1], p[0]))

root = Path(__file__).resolve().parent
header = [
    "id", "group", "radius", "mass", "px", "py", "pz",
    "vx", "vy", "vz", "ox", "oy", "oz", "motion_type",
    "young_modulus", "poisson_ratio", "friction_coeff", "restitution_coeff",
    "energetic", "temperature", "reaction_progress", "specific_heat",
    "thermal_conductivity", "reaction_heat", "arrhenius_prefactor",
    "activation_temperature", "heat_source",
]
with (root / "particles.csv").open("w", newline="", encoding="utf-8") as stream:
    writer = csv.writer(stream)
    writer.writerow(header)
    for idx, (x, y) in enumerate(points, 1):
        writer.writerow([
            idx, 1, 0.45 * SPACING, MASS, x, y, 0,
            0, 0, 0, 0, 0, 0, 1,
            YOUNG_MODULUS, POISSON_RATIO, 0.10, 0.90,
            1, 300.0, 0.0, SPECIFIC_HEAT, THERMAL_CONDUCTIVITY,
            REACTION_HEAT, ARRHENIUS_PREFACTOR, ACTIVATION_TEMPERATURE, 0.0,
        ])

bonds = []
for i, (xi, yi) in enumerate(points):
    for j in range(i + 1, len(points)):
        xj, yj = points[j]
        distance = math.hypot(xj - xi, yj - yi)
        if abs(distance - SPACING) < 1e-9:
            bonds.append((i + 1, j + 1, distance))

with (root / "bonds.csv").open("w", newline="", encoding="utf-8") as stream:
    writer = csv.writer(stream)
    writer.writerow([
        "id", "particle_a_id", "particle_b_id", "rest_length", "kn", "kt",
        "cn", "ct", "active", "conduction_area",
    ])
    for idx, (a, b, distance) in enumerate(bonds, 1):
        writer.writerow([idx, a, b, distance, 1.0e7, 4.0e6, 20.0, 10.0, 1, INTERFACE_AREA])

print(f"generated {len(points)} lattice nodes and {len(bonds)} bonds")
