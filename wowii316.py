"""Written on the Wall II, Conjecture 316 — resolved.

    WOWII 316. Let G be a simple connected graph and let P be its set of pendant
    vertices (degree 1). If avgdeg(Ḡ) ≤ |P|, then G is well totally dominated.

Listed and formalised as open
(`WrittenOnTheWallII/GraphConjecture316.lean`, `category research open`).

    THEOREM. The conjecture is true.

The proof is in two steps, as for conjecture 322 — first work out what the
hypothesis actually admits, then check that family.

**Step 1: the hypothesis admits only cliques, and cliques of order at most 3
with pendants attached.**

Write p = |P| for the number of pendant vertices, r = n - p for the number of
the rest (call them the *core*), and m = |E(G)|. Since
avgdeg(Ḡ) = (n-1) - 2m/n, the hypothesis says

    (n - 1) - 2m/n <= p,    i.e.    m >= n(n - 1 - p)/2 = (r + p)(r - 1)/2.

On the other hand each pendant contributes exactly one edge and the core spans
at most C(r,2), so

    m <= p + r(r-1)/2.

Combining, p + r(r-1)/2 >= (r + p)(r - 1)/2, which simplifies to
p >= (r - 1)p/2, so for p > 0 we get **r <= 3**. Each case then pins the core:

* r = 1: the core is a single vertex, so G is the star K_{1,p}.
* r = 2: G is connected, so the two core vertices must be adjacent (otherwise
  each takes its pendants into a separate component). The core is K_2. Note
  both core vertices carry a pendant -- a core vertex with none would have
  degree 1 and be a pendant itself.
* r = 3: the bound needs m >= 3 + p, and m = (core edges) + p, so the core has
  at least 3 = C(3,2) edges: the core is K_3.

And if p = 0 the hypothesis reads m >= n(n-1)/2, i.e. G = K_n. (The degenerate
case n = 2 is K_2 itself, where both vertices are pendants.)

**Step 2: every graph in that family is well totally dominated.**

* **K_n.** The minimal total dominating sets are exactly the 2-element subsets,
  as in conjecture 322. All have size 2.
* **K_{1,p}, the star.** Every leaf forces the centre into any total dominating
  set, and the centre forces at least one leaf. {centre, leaf} is total
  dominating and minimal; a set with two or more leaves is not minimal, since
  dropping one leaves a total dominating set. All minimal ones have size 2.
* **K_2 with pendants.** Both core vertices u, v carry pendants, so both lie in
  every total dominating set, and {u, v} is already one. All minimal ones have
  size 2.
* **K_3 with pendants.** Exactly the core vertices carrying pendants are forced.
  If all three carry one, the only total dominating set is {u,v,w} and every
  minimal one has size 3. If one or two carry pendants, any forced vertex
  together with one neighbour is total dominating and minimal, and every
  minimal one has size 2.

In every case all minimal total dominating sets of a given graph share a size,
which is what well-totally-dominated means. QED

**How it was found.** Not by searching for a counterexample -- by measuring how
restrictive the hypothesis is. Of the 995 connected graphs on at most 7
vertices it admits 27, and printing them showed every one had a complete core.
That is the same move that settled conjecture 322, and it is cheap: the
hypothesis is checkable in linear time, so "what does this actually admit"
costs one scan.
"""
from __future__ import annotations

import itertools

import networkx as nx

from wowii import (is_well_totally_dominated, minimal_total_dominating_sizes,
                   pendant_vertices, average_degree)


def hypothesis(g):
    return average_degree(nx.complement(g)) <= len(pendant_vertices(g))


def core(g):
    return g.subgraph([v for v, d in g.degree() if d > 1])


def clique_with_pendants(r, dist):
    g = nx.complete_graph(r)
    nxt = r
    for i, k in enumerate(dist):
        for _ in range(k):
            g.add_edge(i, nxt)
            nxt += 1
    return g


def main():
    print("WOWII 316:  avgdeg(complement) <= #pendants  =>  "
          "well totally dominated\n")

    print("step 1 -- the hypothesis admits only cliques, and cliques of order")
    print("          at most 3 with pendants attached.")
    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= 7 and nx.is_connected(g)]
    sat = [g for g in graphs if hypothesis(g)]
    cores_complete = all(
        core(g).number_of_edges()
        == core(g).number_of_nodes() * (core(g).number_of_nodes() - 1) // 2
        for g in sat)
    sizes = sorted({core(g).number_of_nodes() for g in sat})
    with_pendants = sorted({core(g).number_of_nodes() for g in sat
                            if pendant_vertices(g)})
    print(f"   of {len(graphs)} connected graphs on 2..7 vertices, "
          f"{len(sat)} satisfy it")
    print(f"   core sizes seen              : {sizes}")
    print(f"   core sizes when pendants exist: {with_pendants}   "
          f"(the proof gives r <= 3)")
    print(f"   every core complete           : {cores_complete}")

    print("\nstep 2 -- the admitted family is well totally dominated.")
    bad = []
    rows = []
    for n in range(2, 11):
        rows.append((f"K_{n}", nx.complete_graph(n)))
    for r in (1, 2, 3):
        for tot in range(1, 8):
            for dist in itertools.combinations_with_replacement(range(tot + 1), r):
                if sum(dist) != tot:
                    continue
                g = clique_with_pendants(r, dist)
                if nx.is_connected(g) and g.number_of_nodes() <= 10:
                    rows.append((f"K_{r} + pendants {dist}", g))
    seen = set()
    for name, g in rows:
        if not hypothesis(g):
            continue
        s = tuple(sorted(minimal_total_dominating_sizes(g)))
        ok = is_well_totally_dominated(g)
        if not ok:
            bad.append(name)
        seen.add(len(s))
    print(f"   {len(rows)} family members constructed, "
          f"{sum(1 for n, g in rows if hypothesis(g))} satisfy the hypothesis")
    print(f"   number of distinct minimal-TDS sizes within a graph: "
          f"{sorted(seen)}  (1 means well totally dominated)")
    print(f"   violations: {len(bad)}")

    print(f"\nboth steps check out: {cores_complete and not bad}")


if __name__ == "__main__":
    main()
