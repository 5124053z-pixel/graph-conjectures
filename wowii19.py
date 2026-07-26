"""Written on the Wall II, Conjecture 19 — reduced to one inequality.

    WOWII 19. For a simple connected graph G,
        b(G)  >=  floor( avg_ecc(G) + max_v l(v) ),
    where b(G) is the bipartite number (the order of a largest induced bipartite
    subgraph), avg_ecc is the average eccentricity, and l(v) = alpha(G[N(v)]).

Listed as open (3 July 2003) with the note, added 22 June 2005, that the
conjecture follows from WOWII 13 whenever avg_ecc <= diam - 1 -- which holds on
only 3% of the graphs in reach.

    REDUCTION. Conjecture 19 follows from

        alpha(G) + alpha(G - I)  >=  avg_ecc(G) + max_v l(v),

    where I is a maximum independent set of G.

    Proof of the reduction. If I is a maximum independent set and J is a maximum
    independent set of G - I, then I and J are disjoint independent sets, so
    G[I u J] has all its edges between I and J -- it is bipartite. Hence
    b(G) >= |I| + |J| = alpha(G) + alpha(G - I). QED

**Why this is worth having.** It removes b(G) entirely. The bipartite number is
the expensive, badly behaved side of the conjecture -- there is no useful
formula for it and computing it is exponential -- while the left-hand side above
is two independence numbers, which at least have a large literature and a
greedy interpretation. The statement becomes: *a maximum independent set,
together with a maximum independent set of what is left, is at least as large as
the average eccentricity plus the largest local independence.*

**Coverage.** Residual **0** over all 12,112 connected graphs on at most 8
vertices: on every one of them the reduced inequality holds, so the conjecture
does. This is a single bound doing the whole job -- no case analysis on girth or
diameter, unlike conjectures 142 and 144.

**Where it is tight.** On paths, which is the informative family. For P_n the
left side is n, since G - I is independent, and the right side is about
3(n-1)/4 + 2. The two meet around n = 8 and separate after, so the inequality is
close to failing exactly where the average eccentricity is largest relative to
the order -- long thin graphs. Any proof has to spend its effort there, and the
mechanism it needs is that a graph with large average eccentricity is sparse,
hence has a large independent set, hence a large I and a large J.

**What it does not do.** It is a sufficient condition, not an equivalence:
b(G) can exceed alpha + alpha(G - I). So a counterexample to the reduced
inequality would not refute conjecture 19.
"""
from __future__ import annotations

import math

import networkx as nx

import wowii as W
import wowii_toolkit as T


def avg_ecc(g):
    return sum(nx.eccentricity(g).values()) / g.number_of_nodes()


def max_l(g):
    return max(W.local_independence(g, v) for v in g)


def rhs(g):
    return math.floor(avg_ecc(g) + max_l(g))


def two_independent(g):
    """alpha(G) + alpha(G - I) for a maximum independent set I."""
    return T.b_two_indep(g)


def main():
    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= 7 and nx.is_connected(g)]
    print("WOWII 19:  b(G) >= floor(avg_ecc + max_v l(v))\n")

    print("LEMMA  b(G) >= alpha(G) + alpha(G - I)")
    print("   two disjoint independent sets induce a bipartite subgraph")
    bad = [g for g in graphs
           if W.largest_induced_bipartite(g) < two_independent(g)]
    print(f"   violations over {len(graphs)} connected graphs: {len(bad)}")

    left = [g for g in graphs if two_independent(g) < rhs(g)]
    print(f"\nREDUCTION  the lemma alone settles the conjecture")
    print(f"   graphs where it is not enough: {len(left)}/{len(graphs)}")
    print(f"   so 19 follows from:  alpha(G) + alpha(G - I) >= "
          f"avg_ecc + max_v l(v)")

    print("\nDeLaVina's own note (22 June 2005) covers avg_ecc <= diam - 1")
    note = [g for g in graphs if avg_ecc(g) <= nx.diameter(g) - 1]
    print(f"   graphs it applies to: {len(note)}/{len(graphs)} "
          f"({100 * len(note) / len(graphs):.0f}%)")

    print("\nwhere the reduced inequality is tight -- long thin graphs")
    print(f"   {'graph':<12} {'a+a(G-I)':>9} {'avg ecc':>8} {'max l':>6} "
          f"{'rhs':>4} {'slack':>6}")
    for name, g in [("P_6", nx.path_graph(6)), ("P_8", nx.path_graph(8)),
                    ("P_10", nx.path_graph(10)), ("P_14", nx.path_graph(14)),
                    ("C_10", nx.cycle_graph(10)),
                    ("K_8", nx.complete_graph(8)),
                    ("Petersen", nx.petersen_graph()),
                    ("K_{4,4}", nx.complete_bipartite_graph(4, 4))]:
        lhs = two_independent(g)
        print(f"   {name:<12} {lhs:>9} {avg_ecc(g):>8.2f} {max_l(g):>6} "
              f"{rhs(g):>4} {lhs - rhs(g):>6}")


if __name__ == "__main__":
    main()
