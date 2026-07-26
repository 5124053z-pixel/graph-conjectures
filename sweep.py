"""Long-running parameter sweep of the quotient search on bounds 44 and 46.

Prints one line per configuration so progress is visible while it runs, and
keeps the best quotient matrix found for each bound. Anything with negative
slack is printed immediately and in full.
"""
from __future__ import annotations

import argparse
import itertools
import sys
import time

import numpy as np

from quotient_search import (BOUNDS, anneal, realizable, validate, degrees,
                             is_counterexample)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", type=int, nargs="+", default=[2, 3, 4, 5, 6])
    ap.add_argument("--nmax", type=int, nargs="+", default=[40, 150, 600])
    ap.add_argument("--seeds", type=int, default=6)
    ap.add_argument("--iters", type=int, default=40_000)
    ap.add_argument("--restarts", type=int, default=25)
    args = ap.parse_args()

    if validate(verbose=False):
        print("validation against the published counterexamples FAILED")
        raise SystemExit(1)
    print("validated against the twelve published counterexamples", flush=True)

    best = {c: (float("inf"), None, None) for c in BOUNDS}
    t0 = time.time()
    configs = list(itertools.product(sorted(BOUNDS), args.cells, args.nmax,
                                     range(args.seeds)))
    for conj, k, nmax, seed in configs:
        s, n, B = anneal(conj, k, nmax, args.iters, args.restarts,
                         seed=seed * 1013 + k * 17 + nmax, verbose=False)
        tag = ""
        if s < best[conj][0]:
            best[conj] = (s, n, B)
            tag = "  <-- best so far"
        print(f"[{time.time() - t0:7.0f}s] bound {conj}  cells={k}  "
              f"nmax={nmax:>4}  seed={seed}   slack = {s:+.9f}{tag}",
              flush=True)
        if s < -1e-9:
            print(f"\n*** COUNTEREXAMPLE to bound {conj} ***")
            print(f"cell sizes (before scaling): {n}")
            print(f"quotient matrix:\n{B.astype(int)}")
            print(f"degrees {degrees(B)[0]}, neighbour-averages {degrees(B)[1]}")
            print(f"scale factor for realizability: {realizable(n, B)}\n",
                  flush=True)

    print("\n" + "=" * 60)
    for conj, (s, n, B) in sorted(best.items()):
        print(f"bound {conj}: best slack {s:+.9f}")
        if B is not None:
            print(f"  cells {n}\n  quotient\n{B.astype(int)}")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
