"""Written on the Wall II, Conjecture 141 — proved.

    WOWII 141. For a simple connected graph G,
        tree(G)  >=  floor(girth(G)/2) - 1 + max_v l(v),
    where tree(G) is the order of a largest induced tree, l(v) = alpha(G[N(v)])
    is the local independence number, and girth is 0 for an acyclic graph.

Formalised as open in `WrittenOnTheWallII/GraphConjecture141.lean`
(`category research open`).

    THEOREM. The conjecture is true.

Two constructions, one for each range of the girth.

**Case 1: girth <= 5** (including acyclic, where girth is 0 by the convention).
Then floor(g/2) - 1 <= 1, so it is enough to produce an induced tree of order
1 + max_v l(v).

    LEMMA A. tree(G) >= 1 + max_v l(v).
    Proof. Pick v attaining the maximum and let S be a maximum independent set
    of G[N(v)], so |S| = l(v). In G[{v} ∪ S] the only edges join v to S, since S
    is independent, so the induced graph is a star centred at v -- a tree of
    order 1 + l(v). QED

**Case 2: girth g >= 6.** Then G is triangle-free, so N(v) is independent for
every v, hence l(v) = deg(v) and max_v l(v) = Delta.

Put r = floor(g/2) - 1 (so r >= 2), and let v have degree Delta.

    (a) *B(v, r) induces a tree.* Any cycle inside a ball of radius r has length
    at most 2r + 1: a cycle in the induced subgraph contains an edge xy that is
    not in the BFS tree from v, and its fundamental cycle has length at most
    dist(v,x) + dist(v,y) + 1 <= 2r + 1. Here
    2r + 1 = 2*floor(g/2) - 1, which is g - 1 for even g and g - 2 for odd g,
    in both cases strictly less than the girth. So there is no cycle.

    (b) *|B(v, r)| >= Delta + r.* Since B(v,r) induces a tree and G has a cycle,
    B(v,r) is not all of G, so ecc(v) > r and there is a vertex at every
    distance 1, ..., r from v. Counting v, its Delta neighbours and one vertex
    at each distance 2, ..., r gives 1 + Delta + (r - 1) = Delta + r.

Therefore tree(G) >= Delta + r = max_v l(v) + floor(g/2) - 1. QED

**On what is and is not new here.** Both ingredients are standard --
"neighbourhood star" and "balls of small radius are locally tree-like" are
textbook. What appears not to have been done is putting them against this
statement. That is the shape of everything that worked in this project: the
conjecture had been sitting on a list since Graffiti.pc produced it, and nobody
had checked whether two elementary constructions cover its two regimes.

The split point is forced. Lemma A alone gives the conjecture exactly when
floor(g/2) - 1 <= 1, i.e. g <= 5, and the ball construction needs r >= 2, i.e.
g >= 6. The two ranges meet with nothing between them.
"""
from __future__ import annotations

import networkx as nx

from wowii import largest_induced_tree, local_independence


def girth(g):
    return 0 if nx.is_forest(g) else nx.girth(g)


def max_l(g):
    return max(local_independence(g, v) for v in g)


def rhs(g):
    return girth(g) // 2 - 1 + max_l(g)


def ball(g, v, r):
    return list(nx.single_source_shortest_path_length(g, v, cutoff=r).keys())


def star_bound(g):
    """Lemma A's witness: 1 + max_v l(v)."""
    return 1 + max_l(g)


def ball_bound(g):
    """Case 2's witness: the largest ball of radius floor(g/2) - 1 around a
    maximum-degree vertex. Returns None when the girth is below 6."""
    gi = girth(g)
    if gi < 6:
        return None
    r = gi // 2 - 1
    D = max(d for _, d in g.degree())
    return max(len(ball(g, v, r)) for v in g if g.degree(v) == D)


def main():
    print("WOWII 141:  tree(G) >= floor(girth/2) - 1 + max_v l(v)\n")

    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= 7 and nx.is_connected(g)]
    named = [("Petersen", nx.petersen_graph()),
             ("Heawood", nx.heawood_graph()),
             ("Moebius-Kantor", nx.moebius_kantor_graph()),
             ("Pappus", nx.pappus_graph()),
             ("Desargues", nx.desargues_graph()),
             ("C_7", nx.cycle_graph(7)), ("C_9", nx.cycle_graph(9)),
             ("C_12", nx.cycle_graph(12)), ("C_20", nx.cycle_graph(20)),
             ("4x4 grid", nx.grid_2d_graph(4, 4)),
             ("6x6 grid", nx.grid_2d_graph(6, 6)),
             ("Q_4", nx.hypercube_graph(4))]
    named = [(n, nx.convert_node_labels_to_integers(g)) for n, g in named]

    print("LEMMA A  tree(G) >= 1 + max_v l(v)   (the neighbourhood star)")
    bad = [g for g in graphs if largest_induced_tree(g) < star_bound(g)]
    print(f"   violations over {len(graphs)} connected graphs on <= 7 vertices: "
          f"{len(bad)}")
    low = [g for g in graphs if girth(g) <= 5]
    print(f"   settles the conjecture for girth <= 5: {len(low)}/{len(graphs)} "
          f"graphs, all hold: "
          f"{all(largest_induced_tree(g) >= rhs(g) for g in low)}")

    print("\nCASE 2  girth >= 6: the ball of radius floor(g/2) - 1 around a")
    print("        maximum-degree vertex induces a tree of order >= Delta + r")
    print(f"   {'graph':<16} {'g':>2} {'r':>2} {'Delta':>5} {'rhs':>4} "
          f"{'|ball|':>7}  induces a tree?")
    ok = True
    for name, g in named:
        gi = girth(g)
        if gi < 6:
            continue
        r = gi // 2 - 1
        D = max(d for _, d in g.degree())
        trees = all(nx.is_tree(g.subgraph(ball(g, v, r))) for v in g)
        bb = ball_bound(g)
        ok &= trees and bb >= rhs(g)
        print(f"   {name:<16} {gi:>2} {r:>2} {D:>5} {rhs(g):>4} {bb:>7}  {trees}")
    small6 = [g for g in graphs if girth(g) >= 6]
    for g in small6:
        gi = girth(g)
        r = gi // 2 - 1
        trees = all(nx.is_tree(g.subgraph(ball(g, v, r))) for v in g)
        ok &= trees and ball_bound(g) >= rhs(g)
    print(f"   plus the {len(small6)} small graphs of girth >= 6")
    print(f"   every case: ball induces a tree and reaches the bound: {ok}")

    print(f"\nboth cases verified; the two ranges g <= 5 and g >= 6 are "
          f"exhaustive.")


if __name__ == "__main__":
    main()
