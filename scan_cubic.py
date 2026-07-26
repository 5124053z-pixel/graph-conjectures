"""Exhaustive test of  Z(G) <= alpha(G) + 1  over all connected cubic graphs.

The conjecture (TxGraffiti, open since 2017) is that every connected cubic graph
other than K_4 satisfies it. Theorem 18 of arXiv:2410.21724 reduces subcubic
graphs to cubic ones, so this is the whole conjecture.

The scan is only as trustworthy as the enumeration behind it, so
`cubic_graphs.generate` aborts unless its output matches OEIS A002851 at every
level. K_4 must come out as the sole violation; if it does not, something is
wrong with the invariants rather than with the conjecture.
"""
from __future__ import annotations

import argparse
import time

import networkx as nx

from cubic_graphs import generate, KNOWN_COUNTS
from zero_forcing import independence_number, has_forcing_set_of_size, \
    zero_forcing_number


def scan(nmax, verbose=True):
    graphs = generate(nmax, verbose=verbose)
    print()
    violations = []
    tight = {}
    for n in sorted(graphs):
        t0 = time.time()
        n_tight = 0
        for g in graphs[n]:
            a = independence_number(g)
            ok, how = has_forcing_set_of_size(g, a + 1)
            if not ok:
                z = zero_forcing_number(g, lo=a + 2)
                violations.append((n, a, z, sorted(g.edges())))
            else:
                if not has_forcing_set_of_size(g, a, tries=1500)[0]:
                    n_tight += 1
        tight[n] = n_tight
        if verbose:
            print(f"  n = {n:>2}: {len(graphs[n]):>8,} graphs   "
                  f"violations so far: {len(violations)}   "
                  f"Z = alpha+1 exactly: {n_tight:>6,}   "
                  f"[{time.time() - t0:.1f}s]")
    return violations, tight


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=14)
    args = ap.parse_args()

    print("conjecture: connected cubic G != K_4  =>  Z(G) <= alpha(G) + 1")
    print("enumeration verified against OEIS A002851 at every level\n")
    violations, tight = scan(args.n)

    print()
    total = sum(KNOWN_COUNTS[n] for n in KNOWN_COUNTS if n <= args.n)
    print(f"tested every connected cubic graph on <= {args.n} vertices "
          f"({total:,} graphs)")
    if not violations:
        print("no violation at all -- but K_4 must violate it, so the "
              "invariants are suspect")
        return
    for n, a, z, edges in violations:
        tag = " (this is K_4, the known exception)" if n == 4 else \
              "   *** COUNTEREXAMPLE ***"
        print(f"  n={n}  alpha={a}  Z={z}{tag}")
        if n > 4:
            print(f"      edges: {edges}")
    if len(violations) == 1 and violations[0][0] == 4:
        print("\nK_4 is the only violation, as the conjecture asserts.")
        print("No counterexample exists in this range.")


if __name__ == "__main__":
    main()
