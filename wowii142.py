"""Written on the Wall II, Conjecture 142 — proved for girth at most 3.

    WOWII 142. For a simple connected graph G,
        tree(G)  >=  (2/3)*girth(G) + ecc(B),
    where B is the set of vertices of *maximum* eccentricity -- the periphery --
    and ecc(S) = max_v dist(v, S). The girth of a forest is 0.

Listed as open (19 July 2005) and formalised as open in
`WrittenOnTheWallII/GraphConjecture142.lean`.

    THEOREM. The conjecture holds for every connected graph of girth at most 3.

This is the same shape as conjecture 144, and rests on the same kind of fact --
here even simpler, because the periphery is defined by the diameter.

    LEMMA. diam(G) >= 1 + ecc(B) for every connected graph on at least two
    vertices, where B is the periphery.

    Proof. Let v attain k = ecc(B) = dist(v, B).

    If k = 0 the claim is diam >= 1, true for a connected graph on at least
    two vertices.

    If k >= 1 then v is not in B, and B is precisely the set of vertices whose
    eccentricity equals the diameter, so ecc(v) <= diam - 1. Since some vertex
    of B is at distance k from v, k <= ecc(v) <= diam - 1. QED

Compare the corresponding lemma for 144, which needs a two-step argument
because the centre is defined by the *radius*: there one shows ecc(v) > rad
from v not being a centre, and separately k <= rad from the definition of the
radius. For the periphery both halves collapse into one line.

    COROLLARY. If G contains a triangle then

        tree(G) >= diam + 1 >= 2 + ecc(B) = (2/3)*girth + ecc(B),

    the first step because a shortest path is an induced tree. If G is a forest
    then tree(G) = n and the right-hand side is ecc(B) <= n. QED

**What is left.** Girth at least 4, where `(2/3)g` overtakes the 2 that a
shortest path supplies. The gap is small and grows slowly -- at girth 4 the
target exceeds diam + 1 on 34 of the 275 graphs in range, at girth 5 on 21 of
23 -- but it is a genuine gap, and the same one as in 144: every construction
available starts from a shortest cycle, while the conjecture measures its
eccentricity from a distinguished vertex set.
"""
from __future__ import annotations

import networkx as nx

import wowii as W
import wowii_toolkit as T


def girth(g):
    return 0 if nx.is_forest(g) else nx.girth(g)


def periphery_ecc(g):
    e = nx.eccentricity(g)
    hi = max(e.values())
    return W.ecc_set(g, [v for v in e if e[v] == hi], include_S=True)


def rhs(g):
    return (2 / 3) * girth(g) + periphery_ecc(g)


def main():
    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= 7 and nx.is_connected(g)]
    named = [("Petersen", nx.petersen_graph()), ("C_10", nx.cycle_graph(10)),
             ("P_10", nx.path_graph(10)), ("Heawood", nx.heawood_graph()),
             ("K_7", nx.complete_graph(7)), ("star K1,9", nx.star_graph(9)),
             ("5x5 grid", nx.grid_2d_graph(5, 5)),
             ("Q_4", nx.hypercube_graph(4))]
    named = [(n, nx.convert_node_labels_to_integers(g)) for n, g in named]
    everything = graphs + [g for _, g in named]

    print("WOWII 142:  tree(G) >= (2/3)*girth + ecc(B),  B the periphery\n")

    print("LEMMA  diam(G) >= 1 + ecc(B)")
    bad = [g for g in everything if nx.diameter(g) < 1 + periphery_ecc(g)]
    print(f"   violations over {len(everything)} connected graphs: {len(bad)}")
    print(f"   {'graph':<12} {'diam':>5} {'|B|':>4} {'ecc(B)':>7}")
    for name, g in named:
        e = nx.eccentricity(g)
        hi = max(e.values())
        print(f"   {name:<12} {nx.diameter(g):>5} "
              f"{sum(1 for v in e if e[v] == hi):>4} {periphery_ecc(g):>7}")

    print("\nCOROLLARY  the conjecture holds whenever girth <= 3")
    low = [g for g in graphs if girth(g) <= 3]
    print(f"   graphs of girth <= 3 in range: {len(low)}/{len(graphs)}")
    print(f"   all satisfy the conjecture   : "
          f"{all(W.largest_induced_tree(g) >= rhs(g) for g in low)}")
    print(f"   and all fall to diam + 1     : "
          f"{all(nx.diameter(g) + 1 >= rhs(g) for g in low)}")

    hi_g = [g for g in graphs if girth(g) >= 4]
    print(f"\nWHAT IS LEFT  girth >= 4: {len(hi_g)} graphs in range")
    over = [g for g in hi_g if nx.diameter(g) + 1 < rhs(g)]
    print(f"   where diam + 1 is not enough: {len(over)}")
    left = [g for g in hi_g if not any(b(g) >= rhs(g) for b in T.TREE_BOUNDS)]
    print(f"   residual after the whole tree toolkit: {len(left)}")
    print(f"   by girth:")
    for gi in sorted({girth(g) for g in hi_g}):
        sub = [g for g in hi_g if girth(g) == gi]
        o = [g for g in sub if nx.diameter(g) + 1 < rhs(g)]
        print(f"      girth {gi}: {len(sub):>4} graphs, diam+1 insufficient on "
              f"{len(o):>4}")


if __name__ == "__main__":
    main()
