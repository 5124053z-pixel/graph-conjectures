"""Parallel finish of the product-conjecture scan over the 7x7 pairs.

`scan_products.py` walks pairs in order of product size and has already settled
everything with |G|*|H| <= 42 -- 131,279 pairs, i.e. every pair on at most 7
vertices except both being 7 -- with no counterexample. What remains is the
364,231 pairs with both graphs on 7 vertices.

Single-core that is about 17 hours, and the reason is worth noting: the cost is
wildly uneven. A sample of 7x7 pairs runs at 266 pairs/s while the sustained
average is 5.5, because the exact set-cover solve behind gamma and gamma_t is
cheap on sparse products and expensive on dense ones. Uneven cost is exactly the
case where handing chunks to a pool wins, so this splits the work across cores.

Correctness is unchanged -- same `slack`, same exact solvers. Only the schedule
differs.
"""
from __future__ import annotations

import argparse
import itertools
import multiprocessing as mp
import pathlib
import time

import networkx as nx

from domination import (total_domination_number, domination_number,
                        cartesian, direct, validate)
from scan_products import connected_graphs

_GRAPHS: list = []


def _init(graphs):
    global _GRAPHS
    _GRAPHS = graphs


def _slack_of(item):
    k, (i, j) = item
    a, b = _GRAPHS[i], _GRAPHS[j]
    s = total_domination_number(cartesian(a, b)) - domination_number(direct(a, b))
    return k, i, j, s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=max(1, mp.cpu_count() - 2))
    ap.add_argument("--chunk", type=int, default=64)
    ap.add_argument("--checkpoint", default="_prodscan.done",
                    help="file of finished pair indices, so a stopped run resumes")
    args = ap.parse_args()

    if validate(verbose=False):
        print("gamma / gamma_t validation FAILED -- not scanning")
        raise SystemExit(1)
    print("gamma and gamma_t validated")

    by = connected_graphs(7)
    g7 = by[7]
    print(f"{len(g7)} connected graphs on 7 vertices")
    pairs = list(itertools.combinations_with_replacement(range(len(g7)), 2))
    ck = pathlib.Path(args.checkpoint)
    seen = set()
    if ck.exists():
        seen = {int(x) for x in ck.read_text().split()}
        print(f"resuming: {len(seen):,} pairs already done", flush=True)
    todo = [(k, q) for k, q in enumerate(pairs) if k not in seen]
    print(f"{len(pairs):,} pairs total, {len(todo):,} to do, "
          f"{args.workers} workers\n", flush=True)

    t0 = time.time()
    done = tight = 0
    worst = (10 ** 9, None)
    with mp.Pool(args.workers, initializer=_init, initargs=(g7,)) as pool:
        fh = ck.open("a")
        for k, i, j, s in pool.imap_unordered(_slack_of, todo, args.chunk):
            done += 1
            fh.write("%d\n" % k)
            if done % 2000 == 0:
                fh.flush()
            if s == 0:
                tight += 1
            if s < worst[0]:
                worst = (s, (i, j))
            if s < 0:
                print(f"\n*** COUNTEREXAMPLE ***")
                print(f"  G edges: {sorted(g7[i].edges())}")
                print(f"  H edges: {sorted(g7[j].edges())}")
                print(f"  slack = {s}\n", flush=True)
                pool.terminate()
                return
            if done % 5000 == 0:
                el = time.time() - t0
                rate = done / el
                print(f"  [{el:6.0f}s] {done:>7,}/{len(todo):,}  "
                      f"{rate:6.1f} pairs/s  tight {tight:>6,}  "
                      f"min slack {worst[0]:+d}  "
                      f"eta {(len(todo) - done) / rate / 3600:5.2f} h",
                      flush=True)

    print(f"\nfinished {done:,} pairs in {(time.time() - t0) / 3600:.2f} h")
    print(f"no counterexample; minimum slack {worst[0]:+d}")
    print(f"{tight:,} pairs ({100 * tight / done:.1f}%) exactly tight")
    print("\nTogether with the 131,279 pairs already done at smaller product")
    print("sizes, this settles every pair of connected graphs on at most 7")
    print("vertices: 495,510 pairs, no counterexample.")


if __name__ == "__main__":
    mp.freeze_support()
    main()
