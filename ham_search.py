"""Targeted search for a counterexample to the Hamiltonian-path conjectures.

Conjectures 194, 198a, 200 and 217 all conclude "then G has a Hamiltonian path".
A counterexample is therefore a graph with **no** Hamiltonian path that satisfies
the hypothesis, so the search space is not graphs in general -- it is graphs
built to have no Hamiltonian path, of which the cleanest supply is a cut vertex
whose removal leaves three or more components.

Among the 144 connected graphs on 3-7 vertices with no Hamiltonian path, the
best hypothesis slacks are -0.1429 (198a) and -0.2857 (194). Those slacks move
in steps of 1/n, so the granularity improves with n and the gap may close.
"""
from __future__ import annotations

import math
import random

import networkx as nx

import wowii as W


def avg_ecc(g):
    e = nx.eccentricity(g)
    return sum(e.values()) / g.number_of_nodes()


def slack_194(g):
    return (1 + W.avg_l(g)) - W.indep_number(g)


def slack_198(g):
    return (2 + avg_ecc(g)) - W.largest_induced_bipartite(g)


SLACKS = {"194": slack_194, "198a": slack_198}


def blocked(rng, sizes):
    """A cut vertex joined to one vertex of each of several blocks. Removing it
    leaves len(sizes) components, so there is no Hamiltonian path once that is
    at least 3."""
    g = nx.Graph()
    g.add_node("c")
    for i, (k, p) in enumerate(sizes):
        block = [(i, j) for j in range(k)]
        g.add_nodes_from(block)
        for a in range(k):
            for b in range(a + 1, k):
                if rng.random() < p:
                    g.add_edge(block[a], block[b])
        # keep the block connected
        for a in range(1, k):
            if not nx.has_path(g.subgraph(block), block[0], block[a]):
                g.add_edge(block[0], block[a])
        g.add_edge("c", block[0])
    return nx.convert_node_labels_to_integers(g) if nx.is_connected(g) else None


def main():
    rng = random.Random(0)
    best = {k: (-9.0, None) for k in SLACKS}
    for trial in range(200000):
        nblocks = rng.choice([3, 3, 3, 4])
        sizes = [(rng.randint(2, 5), rng.uniform(0.2, 1.0)) for _ in range(nblocks)]
        g = blocked(rng, sizes)
        if g is None or g.number_of_nodes() > 16:
            continue
        if W.has_hamiltonian_path(g):
            continue
        for k, f in SLACKS.items():
            v = f(g)
            if v > best[k][0]:
                best[k] = (v, sorted(g.edges()))
                if v >= 0:
                    print(f"\n*** COUNTEREXAMPLE to #{k}  slack {v:+.4f}")
                    print(f"    edges {sorted(g.edges())}\n", flush=True)
                    return
        if trial % 20000 == 0 and trial:
            print(f"  [{trial:>7,}] " + "  ".join(
                f"#{k} best {best[k][0]:+.4f}" for k in SLACKS), flush=True)
    print("\nfinal: " + "  ".join(f"#{k} best slack {best[k][0]:+.4f}"
                                  for k in SLACKS))
    for k, (v, e) in best.items():
        print(f"  #{k}: {e}")


if __name__ == "__main__":
    main()
