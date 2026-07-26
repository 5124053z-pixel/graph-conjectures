"""Exhaustive test of  gamma_t(G [] H) >= gamma(G x H)  over pairs of small graphs.

The conjecture is tight on every pair tried by hand -- P_3 with P_3, C_4 with
C_5, K_3 with P_4, Petersen with K_2 all give slack exactly 0. A conjecture that
sits on its own boundary everywhere is either a theorem with a clean proof or
false by a small margin, and either way the first thing to do is look at all the
small cases, which nobody appears to have done.

Pairs are enumerated in order of product size so that the cheap evidence arrives
first and the run can be stopped at any point with a complete statement about
everything below some bound.
"""
from __future__ import annotations

import argparse
import itertools
import time

import networkx as nx

from domination import (total_domination_number, domination_number,
                        cartesian, direct, validate)


def connected_graphs(nmax):
    """All connected graphs on 2..nmax vertices, from the atlas (nmax <= 7)."""
    out = {}
    for g in nx.graph_atlas_g():
        n = g.number_of_nodes()
        if 2 <= n <= nmax and nx.is_connected(g):
            out.setdefault(n, []).append(g)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=6, help="largest graph order")
    ap.add_argument("--max-product", type=int, default=10_000)
    args = ap.parse_args()

    if validate(verbose=False):
        print("gamma / gamma_t validation FAILED -- not scanning")
        raise SystemExit(1)
    print("gamma and gamma_t validated against the standard formulas\n")

    by_n = connected_graphs(args.n)
    counts = {n: len(v) for n, v in sorted(by_n.items())}
    print(f"connected graphs by order: {counts}")
    total = sum(counts.values())
    print(f"{total} graphs, {total * (total + 1) // 2:,} unordered pairs\n")

    flat = [(n, i, g) for n in sorted(by_n) for i, g in enumerate(by_n[n])]
    pairs = list(itertools.combinations_with_replacement(flat, 2))
    pairs.sort(key=lambda p: p[0][0] * p[1][0])

    t0 = time.time()
    tested = tight = 0
    worst = (10**9, None)
    for (na, ia, ga), (nb, ib, gb) in pairs:
        if na * nb > args.max_product:
            continue
        gt = total_domination_number(cartesian(ga, gb))
        gd = domination_number(direct(ga, gb))
        s = gt - gd
        tested += 1
        if s == 0:
            tight += 1
        if s < worst[0]:
            worst = (s, (na, ia, nb, ib))
        if s < 0:
            print(f"\n*** COUNTEREXAMPLE ***")
            print(f"  |G| = {na}, |H| = {nb}")
            print(f"  gamma_t(G [] H) = {gt} < gamma(G x H) = {gd}")
            print(f"  G edges: {sorted(ga.edges())}")
            print(f"  H edges: {sorted(gb.edges())}\n", flush=True)
            return
        if tested % 250 == 0:
            print(f"  [{time.time() - t0:6.0f}s] {tested:>6,} pairs   "
                  f"tight (slack 0): {tight:>6,}   min slack {worst[0]:+d}",
                  flush=True)

    print(f"\ntested {tested:,} pairs in {time.time() - t0:.0f}s")
    print(f"no counterexample; minimum slack {worst[0]:+d}")
    print(f"{tight:,} pairs ({100 * tight / max(tested, 1):.1f}%) are exactly tight")


if __name__ == "__main__":
    main()
