"""Cross-entropy-method search for counterexamples to graph-invariant inequalities.

The conjectures targeted here all have the same shape:

    for every (connected) graph G on n vertices,   f(G) >= g(G)

so a counterexample is a single graph with f(G) - g(G) < 0, and the search is
just "minimise the slack". Only the objective changes between conjectures --
the search machinery below is shared. This is the structure Wagner exploits in
"Constructions in combinatorics via neural networks" (arXiv:2104.14516).

A graph on n vertices is encoded as a bit string of length n(n-1)/2, one bit
per possible edge. The cross-entropy method keeps an independent Bernoulli
probability per edge, samples a population, keeps the best-scoring elite
fraction, and moves the probabilities toward the elite's empirical means.

Usage:
    python search.py                 # validate on a conjecture with a known answer
    python search.py --help
"""
from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np
import networkx as nx


# --------------------------------------------------------------------------
# graph encoding
# --------------------------------------------------------------------------

def edge_list(n):
    """The n(n-1)/2 candidate edges, in a fixed order."""
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def build_graph(bits, edges, n):
    g = nx.Graph()
    g.add_nodes_from(range(n))
    g.add_edges_from(e for b, e in zip(bits, edges) if b)
    return g


# --------------------------------------------------------------------------
# invariants
# --------------------------------------------------------------------------

def largest_adjacency_eigenvalue(g):
    a = nx.to_numpy_array(g)
    return float(np.linalg.eigvalsh(a)[-1])


def matching_number(g):
    return len(nx.max_weight_matching(g, maxcardinality=True))


# --------------------------------------------------------------------------
# objectives: return the slack f(G) - g(G); NEGATIVE means counterexample
# --------------------------------------------------------------------------

@dataclass
class Conjecture:
    name: str
    statement: str
    slack: callable          # graph -> float, negative = refuted
    requires_connected: bool = True
    known_answer: str | None = None


def _slack_eigenvalue_plus_matching(g):
    """Conjecture 2.1 in Wagner (from AutoGraphiX, Aouchiche et al.):
    for connected G on n >= 3 vertices,  lambda_1 + mu >= sqrt(n-1) + 1."""
    n = g.number_of_nodes()
    lhs = largest_adjacency_eigenvalue(g) + matching_number(g)
    rhs = math.sqrt(n - 1) + 1.0
    return lhs - rhs


CONJ_2_1 = Conjecture(
    name="wagner-2.1",
    statement="connected G, n>=3:  lambda_1(G) + mu(G) >= sqrt(n-1) + 1",
    slack=_slack_eigenvalue_plus_matching,
    known_answer=(
        "FALSE. Disproved by Stevanovic with a 600-vertex graph; Wagner's "
        "cross-entropy search found a 19-vertex counterexample (a tree). "
        "Used here to validate the search, since the answer is known."
    ),
)

CONJECTURES = {c.name: c for c in (CONJ_2_1,)}


# --------------------------------------------------------------------------
# cross-entropy method
# --------------------------------------------------------------------------

PENALTY = 1e3   # score for graphs that fail the structural precondition


def score_bits(bits, edges, n, conj):
    g = build_graph(bits, edges, n)
    if conj.requires_connected and not nx.is_connected(g):
        return PENALTY, g
    return conj.slack(g), g


def cross_entropy_search(conj, n, iterations=200, population=400,
                         elite_frac=0.08, lr=0.15, seed=0, verbose=True):
    """Minimise conj.slack over graphs on n vertices. Returns (best_score, best_graph)."""
    rng = np.random.default_rng(seed)
    edges = edge_list(n)
    m = len(edges)
    probs = np.full(m, 0.5)
    n_elite = max(2, int(population * elite_frac))

    best_score, best_graph = math.inf, None

    for it in range(1, iterations + 1):
        samples = rng.random((population, m)) < probs
        scored = []
        for row in samples:
            s, g = score_bits(row, edges, n, conj)
            scored.append(s)
            if s < best_score:
                best_score, best_graph = s, g
        scored = np.asarray(scored)

        elite_idx = np.argsort(scored)[:n_elite]
        elite = samples[elite_idx]
        probs = (1 - lr) * probs + lr * elite.mean(axis=0)
        probs = np.clip(probs, 0.01, 0.99)

        if verbose and (it % 20 == 0 or it == 1):
            print(f"  iter {it:4d}   best slack = {best_score:+.6f}"
                  f"   elite mean = {scored[elite_idx].mean():+.4f}")
            if best_score < 0:
                print(f"  --> COUNTEREXAMPLE FOUND at iteration {it}")
                break

    return best_score, best_graph


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--conjecture", default="wagner-2.1", choices=sorted(CONJECTURES))
    ap.add_argument("-n", type=int, default=19, help="number of vertices")
    ap.add_argument("--iterations", type=int, default=200)
    ap.add_argument("--population", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    conj = CONJECTURES[args.conjecture]
    print(f"conjecture : {conj.name}")
    print(f"statement  : {conj.statement}")
    if conj.known_answer:
        print(f"known      : {conj.known_answer}")
    print(f"searching  : n={args.n}, population={args.population}, "
          f"iterations={args.iterations}, seed={args.seed}")
    print()

    score, g = cross_entropy_search(conj, args.n, iterations=args.iterations,
                                    population=args.population, seed=args.seed)

    print()
    print(f"best slack found: {score:+.6f}")
    if score < 0:
        print("RESULT: counterexample found.")
        print(f"  vertices  : {g.number_of_nodes()}")
        print(f"  edges     : {sorted(g.edges())}")
        print(f"  is tree   : {nx.is_tree(g)}")
        print(f"  lambda_1  : {largest_adjacency_eigenvalue(g):.6f}")
        print(f"  mu        : {matching_number(g)}")
        print(f"  threshold : {math.sqrt(args.n - 1) + 1:.6f}")
    else:
        print("RESULT: no counterexample at this n / budget (not a disproof).")


if __name__ == "__main__":
    main()
