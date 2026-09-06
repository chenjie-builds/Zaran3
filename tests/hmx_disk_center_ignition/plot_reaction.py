import csv
from pathlib import Path

import matplotlib.pyplot as plt


iterations = [0, 1000, 2000, 8000]
figure, axes = plt.subplots(1, 4, figsize=(12, 3), constrained_layout=True)
image = None
for axis, iteration in zip(axes, iterations):
    path = Path(f"backup/iter={iteration}/particles.dat")
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    x = [float(row["px"]) * 1000.0 for row in rows]
    y = [float(row["py"]) * 1000.0 for row in rows]
    alpha = [float(row["reaction_progress"]) for row in rows]
    image = axis.scatter(x, y, c=alpha, vmin=0.0, vmax=1.0, cmap="inferno", s=28)
    axis.set_title(f"{iteration * 1.0e-4:.1f} ms")
    axis.set_aspect("equal")
    axis.set_xlabel("x (mm)")
axes[0].set_ylabel("y (mm)")
figure.colorbar(image, ax=axes, label="reaction progress")
figure.savefig("result/reaction_wave.png", dpi=180)
