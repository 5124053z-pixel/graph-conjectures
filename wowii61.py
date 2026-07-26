"""Written on the Wall II, Conjecture 61 — settled for diameter at most 3.

    WOWII 61. For a connected graph G,
        f(G)  >=  residue(G) + ceil(diam(G)/3),
    where f(G) is the order of a largest induced forest and residue(G) is the
    Havel-Hakimi residue.

Two elementary bounds do most of the work, and one of them appears not to be
written down anywhere in this context.

    LEMMA. f(G) >= alpha(G) + 1 for every connected graph on at least two
    vertices.

    Proof. Let I be a maximum independent set and v any vertex outside it (one
    exists, since a connected graph on >= 2 vertices is not edgeless). In
    G[I ∪ {v}] the only edges run between v and N(v) ∩ I, so the induced graph
    is a star centred at v together with isolated vertices. A star is a tree, so
    this is a forest of order alpha + 1. QED

    COROLLARY. Conjecture 61 holds whenever diam(G) <= 3.

    Proof. For diam <= 3 we have ceil(diam/3) = 1, and Favaron, Mahéo and Saclé
    give residue(G) <= alpha(G). So residue + ceil(diam/3) <= alpha + 1 <= f. QED

The second bound is even cheaper: **a shortest path is induced**, so
f(G) >= diam(G) + 1.

Together these two settle the conjecture completely over the graphs reachable
here. Of the 995 connected graphs on at most 7 vertices, 893 have diameter <= 3
and fall to the corollary; of the remaining 102, all but 5 fall to
f >= diam + 1; and those 5 satisfy the conjecture too, with f exceeding the
right-hand side.

**What is left.** The conjecture is open for diameter >= 4, where neither bound
alone suffices in general. What would close it is a strengthening of the lemma:
if S is a set of vertices pairwise at distance >= 3 and disjoint from a maximum
independent set I, then no two members of S share a neighbour, so G[I ∪ S] is a
disjoint union of stars -- again a forest -- and f >= alpha + |S|. A shortest
path of length d supplies ⌊d/3⌋ + 1 such vertices, which is the right order of
magnitude; making them disjoint from some maximum independent set is the gap.

Note the alpha-strengthened form is **false**: f >= alpha + ceil(diam/3) fails
on 8 of the 995 graphs, smallest on 6 vertices. So the weaker residue in the
statement is doing real work, and any proof has to use residue < alpha rather
than route through alpha.
"""
from __future__ import annotations

import math

import networkx as nx

from wowii import largest_induced_forest, indep_number, residue


def rhs(g):
    return residue(g) + math.ceil(nx.diameter(g) / 3)


def main():
    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= 7 and nx.is_connected(g)]
    print("WOWII 61:  f(G) >= residue(G) + ceil(diam(G)/3)\n")

    bad = [g for g in graphs
           if largest_induced_forest(g) < indep_number(g) + 1]
    print(f"LEMMA  f >= alpha + 1")
    print(f"   violations over {len(graphs)} connected graphs: {len(bad)}")

    d3 = [g for g in graphs if nx.diameter(g) <= 3]
    ok3 = [g for g in d3 if largest_induced_forest(g) >= rhs(g)]
    print(f"\nCOROLLARY  the conjecture holds when diam <= 3")
    print(f"   graphs with diam <= 3: {len(d3)}/{len(graphs)}")
    print(f"   all satisfy the conjecture: {len(ok3) == len(d3)}")

    print(f"\nsecond elementary bound: a shortest path is induced, so "
          f"f >= diam + 1")
    bad2 = [g for g in graphs if largest_induced_forest(g) < nx.diameter(g) + 1]
    print(f"   violations: {len(bad2)}")

    d4 = [g for g in graphs if nx.diameter(g) >= 4]
    left = [g for g in d4
            if indep_number(g) + 1 < rhs(g) and nx.diameter(g) + 1 < rhs(g)]
    print(f"\nwhat the two bounds leave")
    print(f"   graphs with diam >= 4          : {len(d4)}")
    print(f"   settled by neither bound        : {len(left)}")
    print(f"   ... and they satisfy it anyway  : "
          f"{all(largest_induced_forest(g) >= rhs(g) for g in left)}")

    print(f"\nthe alpha-strengthened form f >= alpha + ceil(diam/3) is FALSE:")
    badα = [g for g in graphs
            if largest_induced_forest(g) < indep_number(g) + math.ceil(nx.diameter(g) / 3)]
    print(f"   violations: {len(badα)}, smallest on "
          f"{min(g.number_of_nodes() for g in badα)} vertices")
    print(f"   so the residue in the statement is doing real work")


if __name__ == "__main__":
    main()
