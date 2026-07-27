r"""Written on the Wall II, Conjecture 160 — not refuted; proved except for one
narrow case.

    WOWII 160 (DeLaVina's list, 8 August 2005). For a simple connected graph G,
        L_s(G)  >=  max_v l(v)  +  max_v T(v) * c_{C4}(G),
    where L_s is the maximum number of leaves over all spanning trees, l(v) is
    the local independence number, T(v) is the number of triangles containing
    v, and -- **definition 73 of `wowIIdefs.js`, from DeLaVina's own page** --

        c_{C4}(G) is 1 if G is C4-free (NOT necessarily induced), 0 otherwise.

**Two retractions, in order.**

*First*, an earlier version of this file claimed conjecture 160 was false, with
a five-vertex counterexample. That counterexample refutes
`GraphConjecture160.lean`, which reads `c_{C4}` as a **count of induced
4-cycles**. Under that reading the statement fails on 8,985 of the 12,112
connected graphs on at most 8 vertices -- three quarters of them -- which no
conjecture on a twenty-year-old list would do.

*Second*, that file went on to say "no reading of c_C4 rescues it", having
tested a line **labelled** "conjecture 133's convention" that in fact computed
the indicator of an *induced* 4-cycle. Conjecture 133 uses `hasC4`, a 4-cycle as
a subgraph, and the two differ on 310 of the 995 connected graphs on at most 7
vertices.

Definition 73 settles it: the indicator is over 4-cycles **not necessarily
induced**, and under that reading the statement has no violation anywhere it has
been checked. Conjecture 160 is open, and most of it is now proved.

    THEOREM. Conjecture 160 holds for every graph containing a 4-cycle, and for
    every C4-free graph in which some single vertex attains both maxima.

Over the 12,112 connected graphs on at most 8 vertices that leaves **49**.

    CASE 1: G contains a 4-cycle. Then c_{C4} = 0 and the statement is just
    L_s(G) >= max_v l(v). Since l(v) = alpha(G[N(v)]) <= deg(v) <= Delta, and
    L_s >= Delta because a spanning tree rooted at a maximum-degree vertex has
    at least Delta leaves, this is immediate. It covers 11,836 of the 12,112
    graphs in range.

    CASE 2: G is C4-free. Then G[N(v)] is a MATCHING for every v: if x, y, z
    lie in N(v) with x~y and y~z, then x-v-z-y-x is a 4-cycle. Consequently

        T(v) = the number of edges of G[N(v)],
        l(v) = deg(v) - T(v),

    so l(v) + T(v) = deg(v) for every single vertex. If some v attains both
    max_u l(u) and max_w T(w) then the right-hand side is exactly deg(v), which
    is at most Delta <= L_s. That covers 227 of the 276 C4-free graphs in range.

    WHAT IS LEFT: C4-free graphs where the two maxima are attained only at
    different vertices, so the right-hand side is l(u) + T(w) with u != w and
    strictly larger than every individual degree. There are 49 of these on at
    most 8 vertices; all hold, 43 of them with equality.

**The shape of the gap.** With u maximising l and w maximising T, the quantity
l(u) + T(w) exceeds deg(u) = l(u) + T(u) exactly when T(w) > T(u), and exceeds
deg(w) = l(w) + T(w) exactly when l(u) > l(w). So the open case is precisely
"the most locally independent vertex is not the most triangle-heavy one", and
closing it needs a spanning tree that harvests leaves from both neighbourhoods
at once.

**A natural strengthening, and why it fails.** Over all 276 C4-free graphs on at
most 8 vertices, and over Petersen, Heawood, C_9, P_12, balanced trees, stars,
friendship graphs and 26 random C4-free graphs on 9..13 vertices,

    max_v l(v) + max_w T(w)  <=  Delta + 1

holds without exception -- and it would finish the case, since Ls >= Delta + 1
whenever no single vertex attains both maxima. **It is false.** Take a vertex u
of degree 5 with an independent neighbourhood, a vertex w of degree 5 carrying
two disjoint triangles, and one edge joining the two neighbourhoods:

    0-1, 0-2, 0-3, 0-4, 0-5,  6-1, 6-7, 6-8, 6-9, 6-10,  7-8, 9-10

is C4-free on 11 vertices with Delta = 5, max l = l(0) = 5 and max T = T(6) = 2,
so the left side is 7 against Delta + 1 = 6. Conjecture 160 survives it anyway,
with Ls = 8 >= 7.

The moral is the one this project keeps relearning: the strengthening was
supported by every graph in the exhaustive range and by a dozen named families,
and it took ten minutes of *constructing* rather than sampling to break it. A
hill-climbing hunt for a counterexample to conjecture 160 itself, restricted to
C4-free graphs on 9..13 vertices, reaches equality but finds nothing.
"""
from __future__ import annotations

from collections import Counter

import networkx as nx

from wowii import (has_c4_subgraph, local_independence, max_triangles_at_vertex,
                   count_induced_c4, max_leaves_spanning_tree,
                   max_leaves_via_cds)

# K_5 minus the edges 02 and 13 -- the counterexample to the *Lean* statement.
COUNTEREXAMPLE = nx.Graph([(0, 1), (0, 3), (0, 4), (1, 2),
                           (1, 4), (2, 3), (2, 4), (3, 4)])


def max_l(g):
    return max(local_independence(g, v) for v in g)


def readings(g):
    """The ways c_{C4} can be read. Definition 73 picks the last one."""
    n4 = count_induced_c4(g)
    return {
        "count of induced C4  (the Lean file)": n4,
        "indicator, 0 if INDUCED C4 else 1": 0 if n4 else 1,
        "indicator, 1 if induced C4 else 0": 1 if n4 else 0,
        "1 if C4-FREE else 0  (definition 73, the real one)":
            0 if has_c4_subgraph(g) else 1,
    }


def rhs(g):
    """The conjecture's right-hand side, under definition 73."""
    return max_l(g) + max_triangles_at_vertex(g) * (0 if has_c4_subgraph(g) else 1)


def neighbourhood_is_matching(g, v):
    return all(d <= 1 for _, d in g.subgraph(list(g[v])).degree())


def main(nmax=8):
    if nmax <= 7:
        gs = [x for x in nx.graph_atlas_g()
              if 2 <= x.number_of_nodes() <= nmax and nx.is_connected(x)]
    else:
        from allgraphs import connected_graphs
        gs = connected_graphs(nmax, verbose=False)

    print("WOWII 160:  Ls(G) >= max_v l(v) + max_v T(v) * c_C4(G)\n")
    print("definition 73, from DeLaVina's own definitions file:")
    print("   c_C4(G) is 1 if G is C4-free (not necessarily induced), else 0\n")

    print(f"the four readings, over the {len(gs)} connected graphs on 2..{nmax}")
    for name in readings(gs[0]):
        bad = sum(1 for x in gs
                  if max_l(x) + max_triangles_at_vertex(x) * readings(x)[name]
                  > max_leaves_via_cds(x))
        verdict = "HOLDS EVERYWHERE" if not bad else f"{bad} violations"
        print(f"   {name:<52} {verdict}")

    D = lambda g: max(d for _, d in g.degree())
    print("\nCASE 1  G contains a 4-cycle -- the rhs is just max_v l(v)")
    c1 = [g for g in gs if has_c4_subgraph(g)]
    print(f"   {len(c1)} graphs.  l(v) <= deg(v) <= Delta <= Ls, so immediate.")
    print(f"   Ls >= Delta everywhere : "
          f"{all(max_leaves_via_cds(g) >= D(g) for g in gs)}")
    print(f"   the case holds         : "
          f"{all(max_leaves_via_cds(g) >= max_l(g) for g in c1)}")

    print("\nCASE 2  G is C4-free -- every G[N(v)] is a matching")
    c2 = [g for g in gs if not has_c4_subgraph(g)]
    print(f"   {len(c2)} graphs")
    print(f"   every neighbourhood a matching : "
          f"{all(neighbourhood_is_matching(g, v) for g in c2 for v in g)}")
    print(f"   l(v) = deg(v) - T(v)           : "
          f"{all(local_independence(g, v) == g.degree(v) - g.subgraph(list(g[v])).number_of_edges() for g in c2 for v in g)}")
    same = [g for g in c2
            if any(local_independence(g, v) == max_l(g)
                   and g.subgraph(list(g[v])).number_of_edges() == max_triangles_at_vertex(g)
                   for v in g)]
    print(f"   one vertex attains both maxima : {len(same)}   rhs = deg(v) <= Ls")

    rest = [g for g in c2 if g not in same]
    print(f"\nWHAT IS LEFT  {len(rest)} graphs, the two maxima only at "
          f"different vertices")
    print(f"   all hold : {all(max_leaves_via_cds(g) >= rhs(g) for g in rest)}")
    print(f"   slack    : "
          f"{dict(Counter(max_leaves_via_cds(g) - rhs(g) for g in rest))}")

    print("\nthe counterexample to the LEAN statement, for the record")
    g = COUNTEREXAMPLE
    print(f"   K_5 minus 02 and 13:  max l = {max_l(g)}, max T = "
          f"{max_triangles_at_vertex(g)}, induced-C4 count = "
          f"{count_induced_c4(g)}, Ls = {max_leaves_spanning_tree(g)}")
    lean = max_l(g) + max_triangles_at_vertex(g) * count_induced_c4(g)
    print(f"   Lean rhs = {lean} > Ls  -- refutes the file")
    print(f"   real rhs = {rhs(g)} <= Ls  -- the graph has a 4-cycle, so c_C4 = 0")


if __name__ == "__main__":
    main(8)
