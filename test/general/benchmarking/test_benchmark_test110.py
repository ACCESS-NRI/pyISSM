import os
import statistics
import time
from pathlib import Path

import pyissm


# Resolve repo root robustly
THIS_FILE = Path(__file__).resolve()
REPO_ROOT = THIS_FILE.parents[3]  # pyISSM/

ASSETS_DIR = REPO_ROOT / "test" / "assets"


def build_model():
    md = pyissm.model.mesh.triangle(
        pyissm.model.Model(),
        str(ASSETS_DIR / "Exp" / "Square.exp"),
        150000,
    )

    md = pyissm.model.param.set_mask(md, "all", None)

    md = pyissm.model.param.parameterize(
        md,
        str(ASSETS_DIR / "Par" / "SquareShelfConstrained.py"),
    )

    md = pyissm.model.param.set_flow_equation(md, SSA="all")

    md.cluster.np = 3
    md.transient.requested_outputs = ["IceVolume"]

    return md


def run_once():
    md = build_model()
    md = pyissm.model.execute.solve(md, "Transient")
    return md


def test110_runtime_benchmark():
    timings = []
    md = None

    print("Python:", os.sys.executable)
    print("ISSM_DIR:", os.environ.get("ISSM_DIR"))
    print("Repo root:", REPO_ROOT)

    for i in range(3):
        t0 = time.perf_counter()
        md = run_once()
        t1 = time.perf_counter()

        dt = t1 - t0
        timings.append(dt)
        print(f"run {i+1}/3: {dt:.3f} s")

    avg = statistics.mean(timings)
    stdev = statistics.stdev(timings) if len(timings) > 1 else 0.0
    best = min(timings)
    worst = max(timings)

    print("\nPerformance summary for test110:")
    print(f"  avg   : {avg:.3f} s")
    print(f"  stdev : {stdev:.3f} s")
    print(f"  best  : {best:.3f} s")
    print(f"  worst : {worst:.3f} s")

    # sanity check
    assert md is not None
    assert md.results.TransientSolution[2].IceVolume is not None