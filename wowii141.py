"""Written on the Wall II, Conjecture 141 — proved, in its original form.

    WOWII 141 (DeLaVina's list, July 19 2005). For a simple connected graph G,
        tree(G)  >=  (1/2)*girth(G) - 1 + max_v l(v),
    where tree(G) is the order of a largest induced tree and l(v) = alpha(G[N(v)])
    is the local independence number.

**A correction.** An earlier version of this file proved the statement as it
appears in `WrittenOnTheWallII/GraphConjecture141.lean`, which reads

        tree(G) >= floor(girth(G)/2) - 1 + max_v l(v).

For *odd* girth those differ: the original asks for half the girth, the
formalisation for its floor, and the original is stronger by 1 after rounding.
On C_7 the old argument delivered 4 where the original needs 4.5, i.e. 5. So
the old proof established the weaker reading only. Both gaps -- at girth 5 and
at odd girth >= 7 -- are closed below, and the result now covers the original.

    THEOREM. The conjecture is true, in the original (1/2)*girth form.

Write g for the girth, r = floor(g/2) - 1, Delta for the maximum degree. Note
that for g >= 4 the graph is triangle-free, so N(v) is independent and
l(v) = deg(v), giving max_v l(v) = Delta.

    LEMMA A. tree(G) >= 1 + max_v l(v).
    Proof. Pick v attaining the maximum and let S be a maximum independent set
    of G[N(v)]. In G[{v} u S] the only edges join v to S, so it is a star
    centred at v -- a tree of order 1 + l(v). QED

    LEMMA B. If g >= 5 then tree(G) >= 2 + max_v l(v).
    Proof. As above, plus one more step. Girth >= 5 means no two neighbours of
    v share a neighbour other than v, so a vertex w at distance 2 from v is
    adjacent to exactly one member of N(v); adding w to the star keeps it
    induced and a tree. Such a w exists because a graph of girth >= 5 with a
    cycle is not a star. QED

    LEMMA C. If g >= 6 then tree(G) >= Delta + floor(g/2).
    Proof. Take v of degree Delta and consider the ball B(v, r).

    (a) *B(v,r) induces a tree.* A cycle inside a ball of radius r has length
    at most 2r+1: it contains an edge xy outside the BFS tree from v, whose
    fundamental cycle has length at most dist(v,x) + dist(v,y) + 1 <= 2r+1.
    Here 2r+1 = 2*floor(g/2) - 1, which is g-1 for even g and g-2 for odd g --
    below the girth either way.

    (b) *B(v,r) is not all of V*, since it induces a tree while G has a cycle.
    So ecc(v) > r and every level 1, ..., r is nonempty, giving
    |B(v,r)| >= 1 + Delta + (r-1) = Delta + r.

    (c) *One more vertex.* If |B(v,r)| >= Delta + r + 1 we are done, since
    floor(g/2) = r + 1. Otherwise |B(v,r)| = Delta + r, and as levels 2..r
    contribute r-1 vertices in total with each nonempty, **each of those levels
    holds exactly one vertex**. In particular level r is a single vertex z. By
    (b) there is a vertex y at distance r+1 from v; every neighbour of y at
    distance <= r lies at level r, so y is adjacent to z and to nothing else in
    B(v,r). Then B(v,r) u {y} is an induced tree of order Delta + r + 1. QED

The three lemmas cover the three ranges, and they meet exactly:

    girth 0 (acyclic) or 3 or 4:  (1/2)g - 1 <= 1, so Lemma A suffices.
    girth 5:                      (1/2)g - 1 = 1.5, needs 2 + max l -- Lemma B.
    girth >= 6:                   (1/2)g - 1 <= floor(g/2), and max l = Delta,
                                  so Lemma C suffices. QED

**What the correction cost, and why it was invisible.** The floor version and
the original agree for even girth, and graphs of odd girth at least 7 are rare:
there are exactly **two** among the 12,112 connected graphs on at most 8
vertices. An exhaustive check in that range cannot distinguish the two
readings, and the write-up cited the Lean file, so nothing in the loop pointed
at the gap. It surfaced only on reading DeLaVina's original page.
"""
from __future__ import annotations

import networkx as nx

from wowii import largest_induced_tree, local_independence


def girth(g):
    return 0 if nx.is_forest(g) else nx.girth(g)


def max_l(g):
    return max(local_independence(g, v) for v in g)


def rhs_original(g):
    """The statement as DeLaVina published it."""
    return girth(g) / 2 - 1 + max_l(g)


def rhs_lean(g):
    """The statement as formalised, weaker at odd girth."""
    return girth(g) // 2 - 1 + max_l(g)


def ball(g, v, r):
    return list(nx.single_source_shortest_path_length(g, v, cutoff=r))


def bound_star(g):
    """Lemma A."""
    return 1 + max_l(g)


def bound_two_step(g):
    """Lemma B, valid for girth >= 5."""
    return 2 + max_l(g) if girth(g) >= 5 else None


def bound_ball(g):
    """Lemma C, valid for girth >= 6."""
    gi = girth(g)
    if gi < 6:
        return None
    return max(d for _, d in g.degree()) + gi // 2


def main():
    print("WOWII 141:  tree(G) >= (1/2)*girth - 1 + max_v l(v)   (original form)\n")

    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= 7 and nx.is_connected(g)]
    named = [("Petersen", nx.petersen_graph()), ("Heawood", nx.heawood_graph()),
             ("Moebius-Kantor", nx.moebius_kantor_graph()),
             ("Pappus", nx.pappus_graph()), ("Desargues", nx.desargues_graph()),
             ("McGee (girth 7)", nx.LCF_graph(24, [12, 7, -7], 8)),
             ("Tutte-Coxeter (girth 8)",
              nx.LCF_graph(30, [-13, -9, 7, -7, 9, 13], 5)),
             ("C_7", nx.cycle_graph(7)), ("C_9", nx.cycle_graph(9)),
             ("C_13", nx.cycle_graph(13)),
             ("4x4 grid", nx.grid_2d_graph(4, 4)),
             ("Q_4", nx.hypercube_graph(4))]
    named = [(n, nx.convert_node_labels_to_integers(g)) for n, g in named]

    print("the two readings differ only at odd girth")
    odd = [g for g in graphs if girth(g) % 2 == 1 and girth(g) >= 7]
    print(f"   connected graphs on <= 7 vertices with odd girth >= 7: {len(odd)}")
    print(f"   -- which is why an exhaustive check could not see the gap\n")

    print("LEMMA A  tree >= 1 + max l(v)")
    print(f"   violations: {sum(1 for g in graphs if largest_induced_tree(g) < bound_star(g))}")
    print("LEMMA B  girth >= 5  =>  tree >= 2 + max l(v)")
    print(f"   violations: {sum(1 for g in graphs if girth(g) >= 5 and largest_induced_tree(g) < bound_two_step(g))}")
    print("LEMMA C  girth >= 6  =>  tree >= Delta + floor(girth/2)")
    print(f"   violations: {sum(1 for g in graphs if girth(g) >= 6 and largest_induced_tree(g) < bound_ball(g))}")

    print("\nthe structural claims Lemma C rests on")
    bad_tree = bad_level = 0
    for g in graphs + [x for _, x in named]:
        gi = girth(g)
        if gi < 6:
            continue
        r = gi // 2 - 1
        D = max(d for _, d in g.degree())
        for v in g:
            if g.degree(v) != D:
                continue
            lev = {}
            for u, dist in nx.single_source_shortest_path_length(g, v, cutoff=r).items():
                lev.setdefault(dist, []).append(u)
            B = [u for L in lev.values() for u in L]
            if not nx.is_tree(g.subgraph(B)):
                bad_tree += 1
            if len(B) == D + r and any(len(lev.get(i, [])) != 1 for i in range(2, r + 1)):
                bad_level += 1
    print(f"   (a) the ball of radius floor(g/2)-1 induces a tree : {bad_tree == 0}")
    print(f"   (c) when |B| = Delta + r each level 2..r is a single vertex : "
          f"{bad_level == 0}")

    print(f"\n   {'graph':<26} {'g':>2} {'tree':>5} {'original':>9} {'proved':>7}  ok")
    for name, g in named:
        gi = girth(g)
        proved = (bound_ball(g) if gi >= 6
                  else bound_two_step(g) if gi == 5 else bound_star(g))
        print(f"   {name:<26} {gi:>2} {largest_induced_tree(g):>5} "
              f"{rhs_original(g):>9.1f} {proved:>7}  "
              f"{'yes' if proved >= rhs_original(g) else 'NO'}")

    ok = all((bound_ball(g) if girth(g) >= 6
              else bound_two_step(g) if girth(g) == 5
              else bound_star(g)) >= rhs_original(g)
             for g in graphs + [x for _, x in named])
    print(f"\n   the three lemmas cover the original statement everywhere: {ok}")


if __name__ == "__main__":
    main()
