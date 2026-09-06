#!/usr/bin/env python3
"""Validate a permanent inert elastic bond with one fixed endpoint."""

import argparse
import csv
from pathlib import Path


def rows(path: Path) -> dict[int, dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return {int(row["id"]): {key: float(value) for key, value in row.items()}
                for row in csv.DictReader(stream)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sim", type=Path, required=True)
    args = parser.parse_args()
    initial = rows(args.sim / "iter=0" / "particles.dat")
    final = rows(args.sim / "iter=200" / "particles.dat")
    with (args.sim / "iter=200" / "bonds.dat").open(newline="", encoding="utf-8") as stream:
        bond = next(csv.DictReader(stream))

    if abs(final[10]["px"] - initial[10]["px"]) > 1e-12:
        raise AssertionError("fixed bond endpoint moved")
    if final[20]["px"] >= initial[20]["px"] - 1e-4:
        raise AssertionError("stretched bond did not pull the moving endpoint inward")
    if final[20]["vx"] >= 0.0:
        raise AssertionError("moving endpoint has the wrong velocity direction")
    if float(bond["extension"]) <= 0.0 or float(bond["fx"]) <= 0.0:
        raise AssertionError("bond state output is inconsistent with tension")

    print(f"moving endpoint: x={final[20]['px']:.8f}, vx={final[20]['vx']:.8f}")
    print("PASS: permanent inert elastic bond restores its stretched endpoint")


if __name__ == "__main__":
    main()
