import statistics
import time
import pyissm


def build_test110_model():
    md = pyissm.model.mesh.triangle(
        pyissm.model.Model(),
        "../assets/Exp/Square.exp",
        150000,
    )
    md = pyissm.model.param.set_mask(md, "all", None)
    md = pyissm.model.param.parameterize(
        md, "../assets/Par/SquareShelfConstrained.py"
    )
    md = pyissm.model.param.set_flow_equation(md, SSA="all")
    md.cluster.np = 3
    md.transient.requested_outputs = ["IceVolume"]
    return md


def run_test110():
    md = build_test110_model()
    md = pyissm.model.execute.solve(md, "Transient")
    return md


def test110_runtime_benchmark():
    n_runs = 3
    timings = []
    md = None

    for i in range(n_runs):
        t0 = time.perf_counter()
        md = run_test110()
        t1 = time.perf_counter()

        elapsed = t1 - t0
        timings.append(elapsed)
        print(f"run {i + 1}/{n_runs}: {elapsed:.3f} s")

    avg = statistics.mean(timings)
    stdev = statistics.stdev(timings) if len(timings) > 1 else 0.0
    best = min(timings)
    worst = max(timings)

    print("\nPerformance summary for test110:")
    print(f"  avg   : {avg:.3f} s")
    print(f"  stdev : {stdev:.3f} s")
    print(f"  best  : {best:.3f} s")
    print(f"  worst : {worst:.3f} s")

    assert md is not None
    assert md.results.TransientSolution[2].IceVolume is not None