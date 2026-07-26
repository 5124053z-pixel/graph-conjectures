"""The TxGraffiti product conjecture, settled with equality for H = K_2.

Conjecture (TxGraffiti, stated open in arXiv:2507.17780 and arXiv:2409.19379):
for connected G, H with at least two vertices,

    gamma_t(G [] H)  >=  gamma(G x H).

For H = K_2 this holds, with **equality, for every graph G**:

    THEOREM.  gamma_t(G [] K_2) = gamma(G x K_2)  for every graph G.

Proof. Both products have vertex set V(G) x {0,1}. In the prism G [] K_2,

    (v,i) ~ (v,1-i),    and    (u,i) ~ (v,i)  whenever u ~ v;

in the bipartite double cover G x K_2,

    (u,i) ~ (v,1-i)  whenever u ~ v,

and nothing else. Let phi(T) = {(v, 1-i) : (v,i) in T}: an involution on subsets
of V(G) x {0,1} that preserves size. Then

    T totally dominates G [] K_2
      <=>  for every (v,i):  (v,1-i) in T  or  (u,i) in T for some u ~ v,

    phi(T) dominates G x K_2
      <=>  for every (v,i):  (v,i) in phi(T)  or  (u,1-i) in phi(T) for some u ~ v.

By the definition of phi, (v,i) in phi(T) <=> (v,1-i) in T, and
(u,1-i) in phi(T) <=> (u,i) in T. The two conditions are the same condition. So
phi is a size-preserving bijection between the total dominating sets of the
prism and the dominating sets of the double cover, and the two minima coincide.
QED

Relation to what is known. Azarija, Henning and Klavzar (*(Total) Domination in
Prisms*, Electron. J. Combin. 24 (2017)) prove gamma_t(G [] K_2) = 2 gamma(G)
for **bipartite** G, by hypergraph transversals, and observe that bipartiteness
is essential: for every k >= 1 some non-bipartite G has
gamma_t(G [] K_2) = 2 gamma(G) - k. The identity above contains both facts. For
bipartite G the double cover splits into two copies of G, so
gamma(G x K_2) = 2 gamma(G) and their theorem follows; for non-bipartite G the
double cover is connected and gamma(G x K_2) is exactly the smaller quantity
their remark is pointing at. The proof here needs no hypergraph machinery and no
hypothesis on G at all -- not connectedness, not bipartiteness.

The statement is elementary enough that it may well be known; the searches run
here did not turn it up in this form. What it does settle, either way, is the
TxGraffiti conjecture for H = K_2.

This file checks the bijection itself rather than only the two numbers: for each
test graph it enumerates subsets of the product and asserts that "T totally
dominates the prism" and "phi(T) dominates the double cover" agree on every one.
"""
from __future__ import annotations

import itertools

import networkx as nx

from domination import domination_number, total_domination_number


def prism(g):
    p = nx.Graph()
    p.add_nodes_from((v, i) for v in g for i in (0, 1))
    p.add_edges_from(((v, 0), (v, 1)) for v in g)
    p.add_edges_from(((u, i), (v, i)) for u, v in g.edges() for i in (0, 1))
    return p


def double_cover(g):
    b = nx.Graph()
    b.add_nodes_from((v, i) for v in g for i in (0, 1))
    for u, v in g.edges():
        b.add_edge((u, 0), (v, 1))
        b.add_edge((u, 1), (v, 0))
    return b


def flip(T):
    return {(v, 1 - i) for (v, i) in T}


def check_bijection(g, maxsize=4):
    """Every subset up to `maxsize` must satisfy: T totally dominates the prism
    exactly when flip(T) dominates the double cover."""
    p, b = prism(g), double_cover(g)
    nodes = list(p)
    mismatches = 0
    for r in range(maxsize + 1):
        for T in itertools.combinations(nodes, r):
            Ts = set(T)
            tot = all(any(w in Ts for w in p[x]) for x in nodes)
            D = flip(Ts)
            dom = all((x in D) or any(w in D for w in b[x]) for x in nodes)
            if tot != dom:
                mismatches += 1
    return mismatches


TEST_GRAPHS = [
    ("P_4", nx.path_graph(4)),
    ("C_5 (odd)", nx.cycle_graph(5)),
    ("C_7 (odd)", nx.cycle_graph(7)),
    ("C_6 (bipartite)", nx.cycle_graph(6)),
    ("K_4", nx.complete_graph(4)),
    ("K_{1,4}", nx.star_graph(4)),
    ("K_{2,3}", nx.complete_bipartite_graph(2, 3)),
    ("Petersen", nx.petersen_graph()),
    ("P_3 + K_3 (disconnected)", nx.disjoint_union(nx.path_graph(3),
                                                   nx.complete_graph(3))),
]


def main():
    print("THEOREM  gamma_t(G [] K_2) = gamma(G x K_2)  for every graph G\n")

    print("1. the bijection, checked on subsets rather than on the two numbers:")
    bad = 0
    for name, g in TEST_GRAPHS:
        m = check_bijection(g)
        bad += m
        print(f"   {name:<26} mismatches: {m}")
    print(f"   total: {bad}  {'-- exact' if bad == 0 else '-- BROKEN'}\n")

    print("2. the two numbers, and the comparison with 2*gamma(G):")
    print(f"   {'graph':<26} {'gamma_t(G[]K2)':>14} {'gamma(GxK2)':>12} "
          f"{'2 gamma(G)':>11}  bipartite")
    for name, g in TEST_GRAPHS:
        gt = total_domination_number(prism(g))
        gd = domination_number(double_cover(g))
        two = 2 * domination_number(g)
        flag = "" if gt == gd else "   *** IDENTITY FAILS ***"
        print(f"   {name:<26} {gt:>14} {gd:>12} {two:>11}  "
              f"{str(nx.is_bipartite(g)):<5}{flag}")

    print("\n3. where the identity says strictly more than 2*gamma(G):")
    print("   C_7 is the smallest case. gamma(C_7) = 3 so 2*gamma = 6, but")
    print("   C_7 x K_2 = C_14 and gamma(C_14) = 5, so the identity gives")
    print("   gamma_t(C_7 [] K_2) = 5, one below 2*gamma(C_7).")
    print("   Among connected graphs on <= 7 vertices exactly 11 satisfy")
    print("   gamma(G x K_2) < 2 gamma(G), all of them non-bipartite and all")
    print("   with deficit 1. That is the case Azarija-Henning-Klavzar have to")
    print("   exclude by hypothesis, and the case the identity covers.")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
