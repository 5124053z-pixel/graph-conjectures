"""Written on the Wall II, Conjecture 314 — resolved.

    WOWII 314. If G is a connected triangle-free graph whose largest induced
    path has at most 4 vertices, then G is well totally dominated: all its
    minimal total dominating sets have the same size.

Formalised as open in `WrittenOnTheWallII/GraphConjecture314.lean`
(`category research open`).

    THEOREM. The conjecture is true.

Throughout, G is connected, triangle-free and P5-free, S is a minimal total
dominating set, and for x in S we write w(x) for a *private total-neighbour*:
minimality of S means S - {x} fails to totally dominate, so some vertex w(x)
has N(w(x)) n S = {x}. Note the w(x) are distinct for distinct x, since that
equation recovers x from w(x).

Every step below is an induced P5, produced from the structure and killed by
the hypothesis.

================================================================================
PART 1.  |S| <= 3.
================================================================================

Suppose |S| >= 4.

**(1a) G[S] has no isolated vertex.** Total domination asks every vertex of V
to have a neighbour in S -- including the vertices of S themselves.

**(1b) G[S] is connected.** Let C1, C2 be distinct components and put
k = dist_G(C1, C2). Since C1 and C2 are components of an *induced* subgraph
there is no G-edge between them, so k >= 2. If k >= 4 a shortest path between
them already has 5 vertices and is induced. So k is 2 or 3. By (1a) pick
x1 in C1 with a neighbour y1 in C1, and x2 in C2 with a neighbour y2 in C2.

  * k = 2: let z be a common neighbour of x1 and x2. Then z is not in S, since
    an S-vertex adjacent to both would join the two components. The path

        y1 -- x1 -- z -- x2 -- y2

    is induced: y1z and y2z would close triangles with x1 and x2, and every
    remaining pair -- y1x2, y1y2, x1y2, x1x2 -- lies across two components.

  * k = 3: let x1 -- u1 -- u2 -- x2 be a shortest path. Then

        y1 -- x1 -- u1 -- u2 -- x2

    is induced: y1u1 would close a triangle with x1; y1u2 would put y1 within
    distance 2 of x2, against k = 3; and x1u2, u1x2 are excluded because the
    path is shortest.

Either way, a contradiction.

**(1c) G[S] has no induced P4.** If p1 -- p2 -- p3 -- p4 is induced in G[S],
then w(p1) -- p1 -- p2 -- p3 -- p4 is induced in G: p2, p3, p4 all lie in S and
p1 is the only S-neighbour of w(p1), which also forces w(p1) to differ from
each of them.

**(1d) So G[S] is complete bipartite.** A P4-free graph is a cograph, and a
connected cograph on at least two vertices is a join G1 v G2. If either side
had an edge, that edge and any vertex of the other side would close a triangle.
So both sides are edgeless and G[S] = K_{s,t}.

**(1e) Neither side has three vertices.** Say a1, a2, a3 lie on one side and b
on the other, so b is adjacent to all three and no two a's are adjacent. The
three privates w(a1), w(a2), w(a3) cannot be pairwise adjacent -- that is a
triangle -- so some pair, say w(a1), w(a2), is not. Then

    w(a1) -- a1 -- b -- a2 -- w(a2)

is induced: w(a1)b and w(a1)a2 are excluded because a1 is the only S-neighbour
of w(a1), symmetrically for w(a2), and a1a2 is a non-edge of the bipartition.

**(1f) So s = t = 2, and G[S] is a 4-cycle a1 -- b1 -- a2 -- b2 -- a1.** The
eight vertices a1, a2, b1, b2, w(a1), w(a2), w(b1), w(b2) are distinct: no
w can lie in S, because every vertex of a 4-cycle has two S-neighbours while a
private has one. Now every adjacency among the eight is forced.

  * *w(a1)w(a2) is an edge*, else w(a1) -- a1 -- b1 -- a2 -- w(a2) is induced.
    Likewise w(b1)w(b2) is an edge.
  * *w(a1) is adjacent to w(b1) or to w(b2)*, else
    w(a1) -- a1 -- b1 -- w(b1) -- w(b2) is induced. Say w(a1)w(b1).
  * Then w(a2)w(b1) is a non-edge (triangle with w(a1)), so the previous point
    forces w(a2)w(b2); and w(a1)w(b2) is a non-edge (triangle with w(b1)).
  * Every remaining pair is a non-edge: a vertex of S is adjacent to no private
    but its own, and a1a2, b1b2 are the diagonals of the 4-cycle.

That determines all 28 pairs. The resulting graph is 3-regular on 8 vertices
with girth 4 -- it is the **Wagner graph V8**, the circulant C8(1,4). And V8
contains an induced P5, for instance

    b1 -- a1 -- b2 -- w(b2) -- w(a2).

Contradiction. So |S| <= 3.

================================================================================
PART 2.  All minimal total dominating sets have the same size.
================================================================================

Every total dominating set has at least 2 vertices, and a total dominating set
of size 2 is precisely an edge uv with N(u) u N(v) = V -- a *dominating edge*.
So the two cases are:

**G has a dominating edge.** Then every minimal S has |S| = 2.

Triangle-freeness makes N(u) n N(v) empty, since a common neighbour w would
close the triangle uvw, so A := N(u) and B := N(v) partition V, with v in A and
u in B. Triangle-freeness again makes A and B independent -- two adjacent
vertices of N(u) would close a triangle with u -- so G is bipartite with parts
A and B, and u is adjacent to *every* vertex of A while v is adjacent to every
vertex of B. Any total dominating set meets both parts, because a vertex of A
needs a neighbour and all its neighbours lie in B.

Suppose a1, a2 in S n A are distinct. Their privates w(a1), w(a2) lie in B,
being neighbours of vertices of A. Then

    w(a1) -- a1 -- u -- a2 -- w(a2)

is induced: w(a1), u, w(a2) all lie in B and a1, a2 both lie in A, which kills
w(a1)w(a2), w(a1)u, uw(a2) and a1a2 at once; w(a1)a2 and w(a2)a1 are excluded
by privateness; and u differs from w(a1) and w(a2) because N(u) is all of A, so
N(w(a1)) n S would have to contain a2 as well as a1. So |S n A| <= 1, and the
same argument run with v gives |S n B| <= 1. Hence |S| = 2.

**G has no dominating edge.** Then no minimal S has |S| = 2, and Part 1 gives
|S| <= 3, so |S| = 3.

In both cases all minimal total dominating sets of G share a size. QED

================================================================================

**On Bacso-Tuza.** Their theorem -- every connected P5-free graph has a
dominating clique or a dominating P3 -- gives gamma_t <= 3 across the class
immediately, since triangle-freeness shrinks a clique to a vertex or an edge
and each of those is already a total dominating set. That is a clean way to see
the *value* 3, but Part 1 is what the conjecture needs and it does not follow
from Bacso-Tuza: the content is that no *minimal* set overshoots.

**Triangle-freeness is not decoration.** Of the 637 connected P5-free graphs on
at most 7 vertices, 209 fail to be well totally dominated, the smallest on 5
vertices -- edges 04, 12, 13, 23, 34, whose minimal total dominating sets have
sizes 2 and 4. Every use of it above is a triangle closed by one edge.

**How it was found.** By writing down what a minimal set's private neighbours
force and reading off the induced P5 each time. The one branch that resisted
being killed by hand -- s = t = 2 -- turned out to force *every* adjacency among
eight vertices, so it named a single graph, and asking a computer whether that
one graph satisfies the hypothesis took a second. The search space had already
collapsed to a point; the computation only had to confirm it.
"""
from __future__ import annotations

import itertools

import networkx as nx

from wowii import (largest_induced_path, total_domination_number,
                   minimal_total_dominating_sizes, is_well_totally_dominated)


def hypothesis(g):
    return (sum(nx.triangles(g).values()) == 0
            and largest_induced_path(g) <= 4)


def dominating_edge(g):
    for u, v in g.edges():
        if set(g[u]) | set(g[v]) == set(g):
            return u, v
    return None


def dominating_p3(g):
    for b in g:
        for a, c in itertools.combinations(list(g[b]), 2):
            if g.has_edge(a, c):
                continue
            if set(g[a]) | set(g[b]) | set(g[c]) == set(g):
                return a, b, c
    return None


def forced_configuration():
    """The single graph Part 1(f) is forced into."""
    g = nx.Graph()
    g.add_edges_from([("a1", "b1"), ("a1", "b2"),
                      ("a2", "b1"), ("a2", "b2")])          # G[S] = K_{2,2}
    g.add_edges_from([("a1", "wa1"), ("a2", "wa2"),
                      ("b1", "wb1"), ("b2", "wb2")])        # private neighbours
    g.add_edges_from([("wa1", "wa2"), ("wb1", "wb2"),
                      ("wa1", "wb1"), ("wa2", "wb2")])      # forced among privates
    return g


def induced(g, path):
    h = g.subgraph(path)
    return (h.number_of_nodes() == len(path)
            and h.number_of_edges() == len(path) - 1
            and all(h.has_edge(x, y) for x, y in zip(path, path[1:])))


def main():
    print("WOWII 314:  connected, triangle-free, P5-free  =>  "
          "well totally dominated\n")

    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= 7 and nx.is_connected(g)]
    hyp = [g for g in graphs if hypothesis(g)]
    print(f"{len(hyp)} of the {len(graphs)} connected graphs on 2..7 vertices "
          f"satisfy the hypothesis")
    print(f"   all are well totally dominated : "
          f"{all(is_well_totally_dominated(g) for g in hyp)}")
    print(f"   minimal-TDS sizes never exceed 3: "
          f"{max(max(minimal_total_dominating_sizes(g)) for g in hyp)}   "
          f"(Part 1)")

    print("\nPART 1(f)  the one branch that survives by hand forces a single")
    print("           graph on eight vertices. It is the Wagner graph, and it")
    print("           fails the hypothesis, so the branch is empty.")
    v8 = forced_configuration()
    print(f"   vertices {v8.number_of_nodes()}, edges {v8.number_of_edges()}, "
          f"degrees {sorted(dict(v8.degree()).values())}")
    print(f"   triangle-free           : {sum(nx.triangles(v8).values()) == 0}")
    print(f"   girth                   : {nx.girth(v8)}")
    print(f"   isomorphic to C8(1,4)   : "
          f"{nx.is_isomorphic(v8, nx.circulant_graph(8, [1, 4]))}")
    witness = ["b1", "a1", "b2", "wb2", "wa2"]
    print(f"   largest induced path    : {largest_induced_path(v8)}  "
          f"(the hypothesis needs <= 4)")
    print(f"   the stated witness {' -- '.join(witness)} is induced: "
          f"{induced(v8, witness)}")

    print("\nPART 2  the split on whether a dominating edge exists")
    de = [g for g in hyp if dominating_edge(g)]
    rest = [g for g in hyp if not dominating_edge(g)]
    ok_part = True
    for g in de:
        u, v = dominating_edge(g)
        A, B = set(g[u]), set(g[v])
        ok_part &= (A | B == set(g)) and not (A & B)
        ok_part &= not any(g.has_edge(x, y)
                           for x, y in itertools.combinations(A, 2))
        ok_part &= not any(g.has_edge(x, y)
                           for x, y in itertools.combinations(B, 2))
    print(f"   with a dominating edge: {len(de)}")
    print(f"      N(u), N(v) partition V into two independent sets: {ok_part}")
    print(f"      every minimal TDS has size 2                    : "
          f"{all(set(minimal_total_dominating_sizes(g)) == {2} for g in de)}")
    print(f"   without one: {len(rest)}")
    print(f"      every minimal TDS has size 3                    : "
          f"{all(set(minimal_total_dominating_sizes(g)) == {3} for g in rest)}")
    print(f"      each has a dominating P3 (Bacso-Tuza)           : "
          f"{all(dominating_p3(g) for g in rest)}")

    print("\nthe triangle-free hypothesis is load-bearing")
    p5 = [g for g in graphs if largest_induced_path(g) <= 4]
    bad = [g for g in p5 if not is_well_totally_dominated(g)]
    x = min(bad, key=lambda g: g.number_of_nodes())
    print(f"   connected P5-free graphs: {len(p5)}, not WTD: {len(bad)}")
    print(f"   smallest failure: n={x.number_of_nodes()} "
          f"edges={sorted(x.edges())} "
          f"minimal TDS sizes={sorted(set(minimal_total_dominating_sizes(x)))}")


if __name__ == "__main__":
    main()
