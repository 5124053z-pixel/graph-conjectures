"""Local search over cubic graphs, for orders the exhaustive scan cannot reach.

`scan_cubic.py` is definitive but its cost is the number of cubic graphs, which
passes half a million by n = 20. This searches instead, so it reaches n = 40 or
60 -- at the price of never proving absence.

**The objective needs care, and the reason is the mirror image of the Laplacian
failures.** There the quantities were continuous and the trap was that rounding
manufactured counterexamples. Here Z and alpha are integers, so nothing is ever
manufactured -- but `Z - alpha` takes about three distinct values over the whole
search space, giving an objective that is almost everywhere flat. An annealer on
a plateau performs a random walk.

So the objective is layered, each term grading the one above it:

    1000 * (greedy constructions that reach a forcing set of size alpha+1)
    +  10 * (random alpha+1 subsets that force)
    +       (how far random subsets propagate, on average)

The last term is a real number and always moves, so the annealer is never on a
plateau. A state scoring below 1 has defeated both adversaries and is escalated
to the exhaustive test.

Moves preserve 3-regularity: pick two independent edges ab, cd and replace them
with ac, bd, keeping the graph simple and connected.
"""
from __future__ import annotations

import argparse
import math
import random
import time

import networkx as nx

from zero_forcing import (independence_number, forcing_closure,
                          has_forcing_set_of_size, zero_forcing_number, _adj)


def random_cubic(n, rng):
    while True:
        try:
            g = nx.random_regular_graph(3, n, seed=rng.randint(0, 1 << 30))
        except nx.NetworkXError:
            continue
        if nx.is_connected(g):
            return g


def swap(g, rng, tries=60):
    """A 2-opt move that keeps the graph cubic, simple and connected."""
    edges = list(g.edges())
    for _ in range(tries):
        (a, b), (c, d) = rng.sample(edges, 2)
        if len({a, b, c, d}) < 4:
            continue
        if rng.random() < 0.5:
            c, d = d, c
        if g.has_edge(a, c) or g.has_edge(b, d):
            continue
        h = g.copy()
        h.remove_edge(a, b)
        h.remove_edge(c, d)
        h.add_edge(a, c)
        h.add_edge(b, d)
        if nx.is_connected(h):
            return h
    return None


def greedy_hits(adj, nodes, n, k, starts):
    """How many of `starts` greedy constructions reach a forcing set of size k.

    Two revisions got this wrong before it worked.

    The first objective counted only *random* (alpha+1)-subsets that force. It
    reached zero within seconds -- and greedy then found a forcing set anyway,
    every time. Random subsets are a weak adversary, so the search was
    optimising against the wrong opponent and its progress was an illusion.

    The second put greedy in as a penalty applied only when the random count hit
    zero. That built a cliff: states with no random hits were punished by
    +10000, so the annealer learned to avoid exactly the region worth exploring,
    and parked at "one random subset forces" forever. A hard penalty on the
    goal state is a barrier around the goal.

    Counting greedy successes makes it an ordinary graded term instead."""
    hits = 0
    for start in starts:
        s = [start]
        while len(s) < k:
            best, bestsize = None, -1
            for v in nodes:
                if v in s:
                    continue
                size = len(forcing_closure(adj, s + [v]))
                if size > bestsize:
                    best, bestsize = v, size
            s.append(best)
        if len(forcing_closure(adj, s)) == n:
            hits += 1
    return hits


def score(g, samples, rng):
    """Layered objective; smallest is closest to a counterexample. See the
    module docstring for why each layer is there and what broke without it."""
    a = independence_number(g)
    adj, n = _adj(g), g.number_of_nodes()
    nodes = list(g)
    k = a + 1

    hits, reach = 0, 0.0
    for _ in range(samples):
        c = len(forcing_closure(adj, rng.sample(nodes, k)))
        reach += c / n
        if c == n:
            hits += 1
    starts = rng.sample(nodes, min(6, n))
    g_hits = greedy_hits(adj, nodes, n, k, starts)
    return g_hits * 1000 + hits * 10 + reach / samples, a


def search(n, iters=4000, restarts=8, samples=120, seed=0, verbose=True):
    rng = random.Random(seed)
    best = (math.inf, None, None)
    for r in range(restarts):
        g = random_cubic(n, rng)
        cur, a = score(g, samples, rng)
        for it in range(iters):
            T = 6.0 * (1 - it / iters) + 0.05
            h = swap(g, rng)
            if h is None:
                continue
            s, ah = score(h, samples, rng)
            if s < cur or rng.random() < math.exp(-(s - cur) / T):
                g, cur, a = h, s, ah
                if cur < best[0]:
                    best = (cur, g.copy(), a)
                if cur < 1.0:  # no greedy start and no random subset forces
                    ok, how = has_forcing_set_of_size(g, a + 1)
                    if not ok:
                        z = zero_forcing_number(g, lo=a + 2)
                        print(f"\n*** COUNTEREXAMPLE at n = {n} ***")
                        print(f"    alpha = {a}, Z = {z}")
                        print(f"    edges: {sorted(g.edges())}\n", flush=True)
                        return 0, g, a
        if verbose:
            print(f"    restart {r + 1}/{restarts}: best density "
                  f"{best[0]} out of {samples}", flush=True)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, nargs="+", default=[20, 26, 32, 40])
    ap.add_argument("--iters", type=int, default=4000)
    ap.add_argument("--restarts", type=int, default=8)
    ap.add_argument("--samples", type=int, default=120)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    print("searching for cubic G with Z(G) > alpha(G) + 1")
    print("objective: how many random (alpha+1)-subsets force; 0 escalates "
          "to the exhaustive test\n")
    for n in args.n:
        t0 = time.time()
        print(f"n = {n}:")
        d, g, a = search(n, args.iters, args.restarts, args.samples, args.seed)
        if d < 1.0 and g is not None and not has_forcing_set_of_size(g, a + 1)[0]:
            return
        print(f"  best score {d:.4f} at alpha = {a}   "
              f"[{time.time() - t0:.0f}s]   no counterexample\n", flush=True)


if __name__ == "__main__":
    main()
