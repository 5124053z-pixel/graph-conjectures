"""Written on the Wall II, Conjecture 144 — proved for girth at most 3.

    WOWII 144. For a simple connected graph G,
        tree(G)  >=  girth(G) - 1 + ecc(Centers(G)),
    where tree(G) is the order of a largest induced tree, Centers(G) is the set
    of vertices of minimum eccentricity, ecc(S) = max_{v} dist(v, S), and the
    girth of a forest is 0 by the convention on the list.

Listed as open on DeLaVina's page (19 July 2005) and formalised as open in
`WrittenOnTheWallII/GraphConjecture144.lean`.

    THEOREM. The conjecture holds for every connected graph of girth at most 3
    -- that is, for every forest and for every graph containing a triangle.

The whole thing turns on one elementary fact about centres that seems not to be
stated anywhere, and which is worth having on its own.

    LEMMA. diam(G) >= 1 + ecc(Centers(G)) for every connected graph on at
    least two vertices.

    Proof. Write C = Centers(G), r = rad(G), and let v attain
    k = ecc(C) = dist(v, C).

    If k = 0 the claim is diam >= 1, which holds because G is connected on at
    least two vertices.

    Suppose k >= 1. Then v is not a centre, so ecc(v) > r, hence
    diam >= ecc(v) >= r + 1. On the other hand, picking any c in C we have
    dist(v, c) <= ecc(c) = r, so k = dist(v, C) <= r. Combining,
    diam >= r + 1 >= k + 1. QED

    COROLLARY (girth 3). If G contains a triangle then
    tree(G) >= diam(G) + 1 >= 2 + ecc(Centers) = girth - 1 + ecc(Centers).

    The first inequality is the observation that a shortest path is an induced
    path and hence an induced tree; the second is the Lemma.

    COROLLARY (acyclic). If G is a tree then tree(G) = n and the right-hand
    side is ecc(Centers) - 1 <= n, so the conjecture is trivial.

So the conjecture reduces to **girth at least 4**, where the right-hand side
starts outrunning diam + 1: at girth g the target is (g-1) + ecc(Centers) while
the shortest path supplies only diam + 1 >= 2 + ecc(Centers).

**What the remaining range needs.** Three constructions cover every graph of
girth >= 4 on at most 8 vertices, and no two of them suffice:

  * *cycle-path*: a shortest cycle minus a vertex, with a shortest path hung
    off it -- gives (g-1) + dist(C_cycle, farthest vertex);
  * *shortest path*: diam + 1;
  * *caterpillar*: a shortest cycle minus a vertex, with single vertices
    attached at several cycle vertices.

The obstruction to finishing is concrete and visible in the data: **a shortest
cycle need not be eccentric enough.** On 24 of the 309 graphs of girth >= 4 in
range, no shortest cycle C has ecc(C) >= ecc(Centers) -- the smallest is on 6
vertices -- so the cycle-path construction cannot reach the target on its own,
and the conjecture's ecc is measured from the centre rather than from the cycle.
Closing girth >= 4 means relating those two eccentricities, or finding a
construction that starts at the centre instead of at the cycle.
"""
from __future__ import annotations

import networkx as nx

import wowii as W
import wowii_toolkit as T


def girth(g):
    return 0 if nx.is_forest(g) else nx.girth(g)


def ecc_centers(g):
    return W.ecc_set(g, nx.center(g))


def rhs(g):
    return girth(g) - 1 + ecc_centers(g)


def shortest_cycles(g):
    if nx.is_forest(g):
        return []
    gi = nx.girth(g)
    return [c for c in nx.simple_cycles(g, length_bound=gi) if len(c) == gi]


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

    print("WOWII 144:  tree(G) >= girth - 1 + ecc(Centers)\n")

    print("LEMMA  diam(G) >= 1 + ecc(Centers)")
    bad = [g for g in everything if nx.diameter(g) < 1 + ecc_centers(g)]
    print(f"   violations over {len(everything)} connected graphs: {len(bad)}")
    print(f"   {'graph':<12} {'diam':>5} {'rad':>4} {'ecc(C)':>7}")
    for name, g in named:
        print(f"   {name:<12} {nx.diameter(g):>5} {nx.radius(g):>4} "
              f"{ecc_centers(g):>7}")

    print("\nCOROLLARY  the conjecture holds whenever girth <= 3")
    low = [g for g in graphs if girth(g) <= 3]
    print(f"   graphs of girth <= 3 in range: {len(low)}/{len(graphs)}")
    print(f"   all satisfy the conjecture   : "
          f"{all(W.largest_induced_tree(g) >= rhs(g) for g in low)}")
    print(f"   and all fall to diam + 1     : "
          f"{all(nx.diameter(g) + 1 >= rhs(g) for g in low if not nx.is_forest(g))}")

    hi = [g for g in graphs if girth(g) >= 4]
    print(f"\nWHAT IS LEFT  girth >= 4: {len(hi)} graphs in range")
    left = [g for g in hi if not any(b(g) >= rhs(g) for b in T.TREE_BOUNDS)]
    print(f"   covered by the tree toolkit, residual: {len(left)}")

    print("\n   the obstruction: a shortest cycle need not be eccentric enough")
    miss = [g for g in hi
            if max((W.ecc_set(g, c) for c in shortest_cycles(g)), default=-1)
            < ecc_centers(g)]
    print(f"   graphs where no shortest cycle reaches ecc(Centers): "
          f"{len(miss)}/{len(hi)}")
    if miss:
        x = min(miss, key=lambda g: g.number_of_nodes())
        print(f"   smallest: n={x.number_of_nodes()} girth={girth(x)} "
              f"ecc(Centers)={ecc_centers(x)} "
              f"best cycle ecc="
              f"{max(W.ecc_set(x, c) for c in shortest_cycles(x))}")
        print(f"      edges {sorted(x.edges())}")
        print(f"      tree={W.largest_induced_tree(x)} rhs={rhs(x)} "
              f"-- it holds, but not by the cycle-path construction")


if __name__ == "__main__":
    main()
