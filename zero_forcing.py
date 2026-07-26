"""Zero forcing number and independence number, with a validation suite.

Target conjecture (TxGraffiti, open since 2017; see Davila et al.,
arXiv:2410.21724):

    if G is a connected cubic graph and G != K_4, then Z(G) <= alpha(G) + 1.

Theorem 18 of that paper reduces the subcubic case to the cubic one, so cubic
graphs are the whole problem. Proved there for claw-free cubic graphs, for
orientably one-face-embeddable cubic graphs, and asymptotically almost surely
with alpha + 2; an infinite family attains equality. **No exhaustive
computational verification is reported**, which is the gap this file is for.

Two things make this a better target than the spectral conjectures:

* Z and alpha are **integers**. The two worst recurring failures in this
  repository (see the README) were both floating-point equality artifacts on
  conjectures that are exactly tight on an infinite family. This conjecture is
  also exactly tight on an infinite family -- and that now costs nothing,
  because integers have no rounding.
* The automated-refutation literature explicitly *skips* conjectures with
  NP-hard invariants (Roucairol and Cazenave say so directly), so the Monte
  Carlo and deep-RL work that has swept the spectral conjectures has not been
  pointed here.

Definitions. `Z(G)` is the least size of a set S of vertices such that,
starting with S coloured and applying "a coloured vertex with exactly one
uncoloured neighbour colours that neighbour" until no rule applies, every vertex
ends up coloured. `alpha(G)` is the largest independent set.
"""
from __future__ import annotations

import itertools
import random

import networkx as nx


# ---------------------------------------------------------------------------
# zero forcing
# ---------------------------------------------------------------------------

def forcing_closure(adj, seed):
    """Vertices coloured after applying the colour-change rule to exhaustion."""
    coloured = set(seed)
    frontier = list(seed)
    while frontier:
        v = frontier.pop()
        white = [u for u in adj[v] if u not in coloured]
        if len(white) == 1:
            u = white[0]
            coloured.add(u)
            frontier.append(u)
            # v's other neighbours may now have a single white neighbour
            frontier.extend(w for w in adj[u] if w in coloured)
    return coloured


def is_forcing_set(adj, n, seed):
    return len(forcing_closure(adj, seed)) == n


def _adj(g):
    return {v: list(g[v]) for v in g}


def has_forcing_set_of_size(g, k, tries=4000, seed=0, exhaustive=True):
    """Is Z(G) <= k?  Heuristics first, exhaustive search only if they fail.

    The heuristic answers "yes" for the overwhelming majority of graphs in a few
    milliseconds; the exhaustive branch runs only for the rare graph that might
    be a counterexample, which is exactly where the budget should go."""
    adj, n = _adj(g), g.number_of_nodes()
    if k >= n:
        return True, "trivial"
    nodes = list(g)
    rng = random.Random(seed)

    for _ in range(tries):                      # random subsets
        if is_forcing_set(adj, n, rng.sample(nodes, k)):
            return True, "random"

    for start in nodes:                         # greedy: grow from each vertex
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
        if is_forcing_set(adj, n, s):
            return True, "greedy"

    if not exhaustive:
        return False, "heuristics failed"
    for s in itertools.combinations(nodes, k):
        if is_forcing_set(adj, n, s):
            return True, "exhaustive"
    return False, "exhaustive"


def zero_forcing_number(g, lo=None):
    """Exact Z(G). Starts from the minimum degree, which is a valid lower bound."""
    n = g.number_of_nodes()
    lo = lo if lo is not None else max(1, min(dict(g.degree()).values()))
    for k in range(lo, n + 1):
        ok, _ = has_forcing_set_of_size(g, k, tries=600)
        if ok:
            return k
    return n


# ---------------------------------------------------------------------------
# independence number: branch and bound, fast on sparse graphs
# ---------------------------------------------------------------------------

def independence_number(g):
    adj = {v: set(g[v]) for v in g}
    return _mis(adj, set(adj))


def _mis(adj, alive):
    if not alive:
        return 0
    # reduction: a vertex of degree <= 1 is in some maximum independent set
    for v in alive:
        d = adj[v] & alive
        if len(d) <= 1:
            return 1 + _mis(adj, alive - {v} - d)
    v = max(alive, key=lambda x: len(adj[x] & alive))     # branch on max degree
    without = _mis(adj, alive - {v})
    with_v = 1 + _mis(adj, alive - {v} - (adj[v] & alive))
    return max(without, with_v)


# ---------------------------------------------------------------------------
# validation against published values
# ---------------------------------------------------------------------------

def _petersen():
    return nx.petersen_graph()


#: (name, graph, expected Z, expected alpha). Values from the zero-forcing
#: literature; K_4 is the one cubic graph the conjecture excludes, and the
#: Petersen graph and K_{3,3} attain equality, so all three are exactly the
#: cases a wrong implementation would get wrong quietly.
CASES = [
    ("P_5   (path)",        nx.path_graph(5),               1, 3),
    ("C_6   (cycle)",       nx.cycle_graph(6),              2, 3),
    ("K_5   (complete)",    nx.complete_graph(5),           4, 1),
    ("K_4   (excluded!)",   nx.complete_graph(4),           3, 1),
    ("K_{3,3}",             nx.complete_bipartite_graph(3, 3), 4, 3),
    ("K_{2,4}",             nx.complete_bipartite_graph(2, 4), 4, 4),
    ("Petersen",            _petersen(),                    5, 4),
    ("Q_3   (cube)",        nx.convert_node_labels_to_integers(nx.hypercube_graph(3)), 4, 4),
    ("prism K_3 x K_2",     nx.circular_ladder_graph(3),    3, 2),
    ("star K_{1,5}",        nx.star_graph(5),               4, 5),
]


def validate(verbose=True):
    bad = []
    for name, g, z_exp, a_exp in CASES:
        g = nx.convert_node_labels_to_integers(g)
        z, a = zero_forcing_number(g), independence_number(g)
        ok = (z == z_exp) and (a == a_exp)
        if not ok:
            bad.append((name, z, z_exp, a, a_exp))
        if verbose:
            print(f"  {name:<20} Z = {z} (want {z_exp})   "
                  f"alpha = {a} (want {a_exp})   {'OK' if ok else 'MISMATCH'}")
    return bad


def conjecture_holds(g):
    """Z(G) <= alpha(G) + 1.  Returns (holds, Z-or-None, alpha, how)."""
    a = independence_number(g)
    ok, how = has_forcing_set_of_size(g, a + 1)
    return ok, (None if ok else zero_forcing_number(g, lo=a + 2)), a, how


if __name__ == "__main__":
    print("validating Z and alpha against published values\n")
    bad = validate()
    print()
    if bad:
        print(f"MISMATCHES: {bad}")
        raise SystemExit(1)
    print("all published values reproduced\n")
    print("the conjecture on the validation cases (cubic ones only):")
    for name, g, _, _ in CASES:
        g = nx.convert_node_labels_to_integers(g)
        degs = set(dict(g.degree()).values())
        if degs != {3}:
            continue
        ok, z, a, how = conjecture_holds(g)
        mark = "holds" if ok else f"FAILS  (Z = {z} > alpha + 1 = {a + 1})"
        print(f"  {name:<20} alpha = {a}   {mark}")
