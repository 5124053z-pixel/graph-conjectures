"""Written on the Wall II, Conjecture 100 — reduced to a bounded range.

    WOWII 100. For a connected graph G whose complement is connected,
        alpha(G)  <=  ceil( ( max_v l(v) + 0.5 * ||deg(Gc)||_2 ) / 2 ),
    where l(v) = alpha(G[N(v)]) is the local independence number, Gc is the
    complement, and ||deg(Gc)||_2 = sqrt( sum_v deg_Gc(v)^2 ).

Formalised as open in `WrittenOnTheWallII/GraphConjecture100.lean`
(`category research open`).

    THEOREM. The conjecture holds for every graph with alpha(G) >= 16.

So it can only fail in the range alpha <= 15.

Two one-line lemmas, and the arithmetic closes.

    LEMMA 1. ||deg(Gc)||_2 >= (alpha - 1) * sqrt(alpha).

    Proof. A maximum independent set I of G is a *clique* of size alpha in Gc,
    so every one of its alpha vertices has Gc-degree at least alpha - 1. Hence
    sum_v deg_Gc(v)^2 >= alpha * (alpha - 1)^2. QED

    LEMMA 2. max_v l(v) >= 2 for every connected non-complete graph.

    Proof. G is not complete, so some pair is non-adjacent; a shortest path
    between such a pair has at least three vertices, and its first three give
    x, z, y with x, y both in N(z) and xy a non-edge. So l(z) >= 2. QED

Combining, and using ceil(t) >= t,

    RHS  >=  ( 2 + 0.5 * (alpha - 1) * sqrt(alpha) ) / 2
         =   1 + (alpha - 1) * sqrt(alpha) / 4.

The conjecture therefore holds as soon as

    alpha <= 1 + (alpha - 1) * sqrt(alpha) / 4
      <=>  (alpha - 1) <= (alpha - 1) * sqrt(alpha) / 4
      <=>  4 <= sqrt(alpha)
      <=>  alpha >= 16.

(Note alpha >= 16 covers alpha >= 2, so Lemma 2 applies.) QED

**Why this conjecture resisted everything else.** It is the only one of the
nineteen that bounds a hard invariant from *above*. Every other technique in
this project works by exhibiting a witness -- an induced star, a ball of small
radius, a greedy induced path, a connected dominating set -- and a witness
proves a maximum-type invariant is **large**, never that it is small. Bounding
alpha from above means ruling out *all* independent sets, which is computing
alpha. Four valid free upper bounds were tried (Gallai's n - nu, n - delta,
the clique bound (1+sqrt(1+8*mbar))/2, and n/2 under a perfect matching) and
together they left 98 of the 739 applicable graphs on at most 7 vertices
untouched, against 100 for Gallai alone. The other three add almost nothing.

What works instead is not a bound on alpha but an argument that makes the
*right-hand side* grow with alpha -- which it does, because a large independent
set is a large clique in the complement and therefore forces complement degrees
up. That is the same move as everywhere else in this project, applied to the
other side of the inequality.

**What the remaining range looks like.** For alpha <= 15 there is a real
trade-off the proof above throws away. Write t for the largest number of
I-neighbours of a vertex outside I. Then max_v l(v) >= t, since N(u) n I is
independent; and each vertex outside I has Gc-degree at least alpha - t, so
||deg(Gc)||_2^2 >= alpha*(alpha-1)^2 + (n-alpha)*(alpha-t)^2. Large t helps the
first term of the right-hand side, small t helps the second. Exploiting that
seesaw is what would push the threshold below 16.
"""
from __future__ import annotations

import math

import networkx as nx

from wowii import indep_number, local_independence, degree_l2_norm, conj_100


def rhs(g):
    return math.ceil((max(local_independence(g, v) for v in g)
                      + 0.5 * degree_l2_norm(nx.complement(g))) / 2)


def main():
    print("WOWII 100:  alpha <= ceil((max l(v) + 0.5*||deg(Gc)||_2) / 2)\n")

    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= 7 and nx.is_connected(g)]

    print("LEMMA 1  ||deg(Gc)||_2 >= (alpha - 1) * sqrt(alpha)")
    bad = [g for g in graphs
           if degree_l2_norm(nx.complement(g))
           < (indep_number(g) - 1) * math.sqrt(indep_number(g)) - 1e-9]
    print(f"   violations over {len(graphs)} connected graphs: {len(bad)}")

    print("\nLEMMA 2  max_v l(v) >= 2 for a connected non-complete graph")
    nc = [g for g in graphs if indep_number(g) >= 2]
    bad2 = [g for g in nc if max(local_independence(g, v) for v in g) < 2]
    print(f"   non-complete connected graphs: {len(nc)}, violations: {len(bad2)}")

    print("\nTHRESHOLD  alpha <= 1 + (alpha-1)*sqrt(alpha)/4  <=>  alpha >= 16")
    print(f"   {'alpha':>6} {'1 + (a-1)sqrt(a)/4':>20}  holds?")
    for a in (2, 4, 8, 12, 14, 15, 16, 17, 20):
        v = 1 + (a - 1) * math.sqrt(a) / 4
        print(f"   {a:>6} {v:>20.3f}  {'yes' if a <= v else 'no'}")
    first = next(a for a in range(2, 100) if a <= 1 + (a - 1) * math.sqrt(a) / 4)
    print(f"   first alpha at which it holds: {first}")

    print("\nthe theorem says nothing in reach -- alpha >= 16 needs n >= 16 --")
    print("so the conjecture is also checked directly where it can be:")
    ap = [g for g in graphs if conj_100(g)[0]]
    print(f"   applicable on 2..7 vertices: {len(ap)}, "
          f"all hold: {all(conj_100(g)[1] for g in ap)}")
    tight = [g for g in ap if indep_number(g) == rhs(g)]
    print(f"   tight (alpha exactly equal to the bound): {len(tight)}")
    print(f"   maximum alpha seen in range: {max(indep_number(g) for g in ap)} "
          f"-- far below the threshold, as expected")

    print("\nthe seesaw the proof throws away")
    print(f"   {'graph':<28} {'alpha':>5} {'t':>3} {'maxl':>5} {'||degGc||':>10}")
    for name, g in [("Petersen", nx.petersen_graph()),
                    ("K_{5,5}", nx.complete_bipartite_graph(5, 5)),
                    ("C_16", nx.cycle_graph(16)),
                    ("P_16", nx.path_graph(16)),
                    ("Q_4", nx.convert_node_labels_to_integers(
                        nx.hypercube_graph(4)))]:
        a = indep_number(g)
        I = max(nx.find_cliques(nx.complement(g)), key=len)
        t = max((len(set(g[u]) & set(I)) for u in g if u not in I), default=0)
        print(f"   {name:<28} {a:>5} {t:>3} "
              f"{max(local_independence(g, v) for v in g):>5} "
              f"{degree_l2_norm(nx.complement(g)):>10.2f}")


if __name__ == "__main__":
    main()
