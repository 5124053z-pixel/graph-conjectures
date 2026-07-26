"""Exhaustive generation of all connected graphs, count-verified.

networkx's atlas stops at 7 vertices, which has been the ceiling on every
exhaustive scan in this repository. This lifts it.

The rule is the obvious one: every graph on n vertices is a graph on n-1
vertices plus a new vertex joined to some subset of the old ones. Applied to a
complete list on n-1 vertices it produces a complete list on n, and isomorphic
duplicates are rejected.

As with `cubic_graphs.py`, the rule is checked rather than trusted: the number
of connected graphs produced must match OEIS A001349 at every level. A silently
incomplete enumeration would turn "no counterexample" into a false statement and
is undetectable from inside the scan, so the run aborts on the first
disagreement.

Cost is dominated by isomorphism rejection, so candidates are bucketed by a
cheap invariant first -- degree sequence together with the sorted multiset of
closed-walk counts -- and `is_isomorphic` runs only inside a bucket.
"""
from __future__ import annotations

import argparse
import itertools
import time

import numpy as np
import networkx as nx

#: OEIS A001349, connected graphs on n unlabelled vertices.
CONNECTED_COUNTS = {1: 1, 2: 1, 3: 2, 4: 6, 5: 21, 6: 112, 7: 853,
                    8: 11117, 9: 261080}
#: OEIS A000088, all graphs on n unlabelled vertices.
ALL_COUNTS = {1: 1, 2: 2, 3: 4, 4: 11, 5: 34, 6: 156, 7: 1044,
              8: 12346, 9: 274668}


def certificate(g, n):
    """Cheap isomorphism-invariant bucket key: degree sequence plus closed-walk
    counts. Only a bucketing device -- an exact test runs inside each bucket, so
    a weak key costs time and never costs correctness."""
    a = nx.to_numpy_array(g, nodelist=sorted(g)).astype(np.int64)
    a2 = a @ a
    a3 = a2 @ a
    return (tuple(sorted(int(x) for x in a.sum(axis=1))),
            tuple(sorted(int(x) for x in np.diag(a2))),
            tuple(sorted(int(x) for x in np.diag(a3))),
            int(a3.sum()))


def extend(g, n):
    """All graphs on n vertices obtained by adding one vertex to g."""
    old = sorted(g)
    for r in range(len(old) + 1):
        for nbrs in itertools.combinations(old, r):
            h = g.copy()
            h.add_node(n - 1)
            h.add_edges_from((n - 1, u) for u in nbrs)
            yield h


def generate(nmax, verbose=True):
    """{n: [all graphs on n vertices]}, complete and up to isomorphism."""
    level = [nx.empty_graph(1)]
    out = {1: level}
    for n in range(2, nmax + 1):
        t0 = time.time()
        buckets = {}
        for g in level:
            for h in extend(g, n):
                key = certificate(h, n)
                bucket = buckets.setdefault(key, [])
                if not any(nx.is_isomorphic(h, k) for k in bucket):
                    bucket.append(h)
        level = [h for b in buckets.values() for h in b]
        out[n] = level
        conn = sum(1 for h in level if nx.is_connected(h))
        exp_all, exp_conn = ALL_COUNTS.get(n), CONNECTED_COUNTS.get(n)
        ok = ((exp_all is None or len(level) == exp_all)
              and (exp_conn is None or conn == exp_conn))
        if verbose:
            print(f"  n = {n}: {len(level):>7,} graphs "
                  f"(expect {exp_all if exp_all else '?':>7}), "
                  f"{conn:>7,} connected (expect {exp_conn if exp_conn else '?':>7})"
                  f"   {'OK' if ok else '*** MISMATCH ***'}   "
                  f"[{time.time() - t0:.1f}s]")
        if not ok:
            raise RuntimeError(
                f"generation is incomplete at n={n}. The scans must not be "
                f"trusted until this is fixed.")
    return out


def connected_graphs(nmax, verbose=False):
    """All connected graphs on 2..nmax vertices, count-verified."""
    lv = generate(nmax, verbose=verbose)
    return [g for n in range(2, nmax + 1) for g in lv[n] if nx.is_connected(g)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=8)
    args = ap.parse_args()
    print("generating all graphs, verifying counts against OEIS A000088 "
          "and A001349\n")
    generate(args.n)
    print("\nall counts match; the enumeration is complete in this range")


if __name__ == "__main__":
    main()
