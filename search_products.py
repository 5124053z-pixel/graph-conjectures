"""Local search over pairs (G, H) for a counterexample to the product conjecture.

`scan_products.py` is exhaustive but only reaches both graphs on <= 7 vertices.
This searches over pairs at fixed orders, so it can look at asymmetric shapes --
a small G against a larger H and vice versa -- that the exhaustive scan never
sees.

The objective is the slack itself, `gamma_t(G [] H) - gamma(G x H)`. Unlike the
cubic search this needs no layering: slack 0 is common (7.5% of pairs on <= 6
vertices) but slack 1, 2, 3 are all well populated too, so the landscape has
real structure rather than one enormous plateau. Moves that leave the slack
unchanged are accepted freely, which turns the slack-0 set into exploration
rather than a trap -- the lesson from failure 6 and failure 10, applied in
advance for once.

Cost is dominated by exact domination on the products, which grows sharply:
a 40-vertex product takes tenths of a second, a 60-vertex one can take two
minutes. Product size is capped accordingly.
"""
from __future__ import annotations

import argparse
import math
import random
import time

import networkx as nx

from domination import (total_domination_number, domination_number,
                        cartesian, direct, validate)


def random_connected(n, rng, p=0.35):
    while True:
        g = nx.gnp_random_graph(n, p, seed=rng.randint(0, 1 << 30))
        if n >= 2 and nx.is_connected(g):
            return g


def perturb(g, rng, tries=40):
    """Flip one edge, keeping the graph connected and non-empty."""
    n = g.number_of_nodes()
    for _ in range(tries):
        u, v = rng.sample(range(n), 2)
        h = g.copy()
        if h.has_edge(u, v):
            h.remove_edge(u, v)
        else:
            h.add_edge(u, v)
        if h.number_of_edges() and nx.is_connected(h):
            return h
    return None


def slack(g, h):
    return total_domination_number(cartesian(g, h)) - domination_number(direct(g, h))


def search(ng, nh, iters, seed, verbose=True):
    rng = random.Random(seed)
    g, h = random_connected(ng, rng), random_connected(nh, rng)
    cur = slack(g, h)
    best = cur
    for it in range(iters):
        which = rng.random() < 0.5
        cand = perturb(g if which else h, rng)
        if cand is None:
            continue
        g2, h2 = (cand, h) if which else (g, cand)
        s = slack(g2, h2)
        T = 1.5 * (1 - it / iters) + 0.02
        if s <= cur or rng.random() < math.exp(-(s - cur) / T):
            g, h, cur = g2, h2, s
            if cur < best:
                best = cur
            if cur < 0:
                print(f"\n*** COUNTEREXAMPLE  |G|={ng} |H|={nh}  slack={cur}")
                print(f"    G: {sorted(g.edges())}")
                print(f"    H: {sorted(h.edges())}\n", flush=True)
                return cur, g, h
    return best, g, h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--iters", type=int, default=250)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--max-product", type=int, default=44)
    args = ap.parse_args()

    if validate(verbose=False):
        print("gamma / gamma_t validation FAILED -- not searching")
        raise SystemExit(1)
    print("gamma and gamma_t validated\n")

    shapes = [(a, b) for a in range(3, 12) for b in range(3, 12)
              if a <= b and a * b <= args.max_product]
    print(f"{len(shapes)} shapes, {args.seeds} seeds each, "
          f"{args.iters} iterations\n")
    for ng, nh in shapes:
        t0 = time.time()
        best = 99
        for seed in range(args.seeds):
            s, g, h = search(ng, nh, args.iters, seed)
            if s < 0:
                return
            best = min(best, s)
        print(f"  |G|={ng:>2} |H|={nh:>2} (product {ng * nh:>3})   "
              f"min slack {best:+d}   [{time.time() - t0:.0f}s]", flush=True)


if __name__ == "__main__":
    main()
