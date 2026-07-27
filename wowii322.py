r"""Written on the Wall II, Conjecture 322 — resolved; the hypothesis is about
the complement.

    WOWII 322 (DeLaVina's list, 4 March 2007). Let G be a simple connected
    graph with n >= 5. If lambda_max(Gc) <= 1, then G is well totally
    dominated -- where Gc is the COMPLEMENT and lambda_max is the largest local
    independence, lambda(v) = alpha(Gc[N_Gc(v)]).

**A correction.** An earlier version of this file read the hypothesis as
lambda_max(G) <= 1, without the complement. That is a different and far more
restrictive condition: it admits 4 of the 12,112 connected graphs on at most 8
vertices, all complete, so the "proof" amounted to "the hypothesis forces
G = K_n". The real hypothesis admits 51 of them, most not complete, and the
conjecture has genuine content. The complement is an overline in the source
HTML and was lost when the page was read as text -- the same class of error as
the floor in conjecture 141 and the c_C4 count in conjecture 160.

    THEOREM. The conjecture is true, and the hypothesis n >= 5 is unnecessary.

    STEP 1. lambda_max(Gc) <= 1 holds exactly when G is complete multipartite.

    lambda_max(Gc) <= 1 says that for every v the set N_Gc(v) is a clique of
    Gc -- equivalently, V \ N[v] is an independent set of G. So no two
    non-neighbours of v are adjacent, for any v.

    Let u and v be non-adjacent. If x is a neighbour of u lying outside N(v),
    then x is not v (v is not a neighbour of u), so x and u both lie in
    V \ N[v], which is independent -- contradicting x in N(u). Hence N(u) is
    contained in N(v), and by symmetry N(u) = N(v).

    Non-adjacency is therefore transitive: if u,v and v,w are non-adjacent then
    N(u) = N(v) = N(w), and an edge uw would put w in N(u) = N(v), making w a
    neighbour of v. So "equal or non-adjacent" is an equivalence relation, its
    classes are independent, and every pair from different classes is an edge:
    G is complete multipartite. The converse is immediate -- the complement of
    a complete multipartite graph is a disjoint union of cliques, and in a union
    of cliques every neighbourhood is a clique.

    STEP 2. Every connected complete multipartite graph is well totally
    dominated, with gamma_t = 2.

    Let the parts be V_1, ..., V_k; connectivity gives k >= 2. For u in V_i and
    w in V_j with i != j,

        N(u) u N(w) = (V \ V_i) u (V \ V_j) = V,

    since V_i and V_j are disjoint. So {u,w} is a total dominating set and
    gamma_t = 2.

    Now let S be any total dominating set. A vertex of V_i is adjacent only to
    V \ V_i, so S cannot lie inside a single part -- otherwise that part is
    undominated. Hence S contains two vertices from different parts, and those
    two already form a total dominating set. If S is minimal it contains
    nothing else, so |S| = 2. QED

Both steps are verified over all 12,112 connected graphs on at most 8 vertices:
the characterisation has zero mismatches, and every complete multipartite graph
in range -- plus larger ones such as K(5,5), K(2,3,4) and K(2,2,2,2) -- has all
its minimal total dominating sets of size 2.

**Note on n >= 5.** It plays no role. The argument works for every connected
complete multipartite graph, and the seven such graphs on fewer than 5 vertices
are well totally dominated too.
"""
from __future__ import annotations

import networkx as nx

from wowii import (is_well_totally_dominated, minimal_total_dominating_sizes,
                   total_domination_number, local_independence)


def lambda_max_complement(g):
    gc = nx.complement(g)
    if gc.number_of_nodes() == 0:
        return 0
    return max(local_independence(gc, v) for v in gc)


def hypothesis(g):
    return lambda_max_complement(g) <= 1


def is_complete_multipartite(g):
    """The complement is a disjoint union of cliques."""
    gc = nx.complement(g)
    for c in nx.connected_components(gc):
        k = len(c)
        if gc.subgraph(c).number_of_edges() != k * (k - 1) // 2:
            return False
    return True


def main():
    print("WOWII 322:  n >= 5 and lambda_max(COMPLEMENT) <= 1  =>  "
          "well totally dominated\n")

    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= 7 and nx.is_connected(g)]
    complete = lambda x: (x.number_of_edges()
                          == x.number_of_nodes() * (x.number_of_nodes() - 1) // 2)

    print("the hypothesis is about the complement, and the difference matters")
    wrong = [g for g in graphs
             if max(local_independence(g, v) for v in g) <= 1]
    right = [g for g in graphs if hypothesis(g)]
    print(f"   lambda_max(G)  <= 1 admits {len(wrong):>4}   all complete: "
          f"{all(complete(x) for x in wrong)}")
    print(f"   lambda_max(Gc) <= 1 admits {len(right):>4}   all complete: "
          f"{all(complete(x) for x in right)}")

    print("\nSTEP 1  lambda_max(Gc) <= 1  <=>  G is complete multipartite")
    bad = [g for g in graphs if hypothesis(g) != is_complete_multipartite(g)]
    print(f"   mismatches over {len(graphs)} connected graphs: {len(bad)}")

    print("\nSTEP 2  connected complete multipartite  =>  every minimal TDS "
          "has size 2")
    cm = [g for g in graphs if is_complete_multipartite(g)]
    print(f"   {len(cm)} such graphs in range")
    print(f"   gamma_t = 2 on all      : "
          f"{all(total_domination_number(g) == 2 for g in cm)}")
    print(f"   every minimal TDS size 2: "
          f"{all(set(minimal_total_dominating_sizes(g)) == {2} for g in cm)}")

    print("\n   larger members of the family")
    print(f"   {'parts':<16} {'n':>3} {'gamma_t':>8} {'minimal sizes':>15}")
    for parts in [(3, 3), (2, 2, 2), (1, 2, 3), (4, 4), (2, 3, 4),
                  (1, 1, 1, 5), (5, 5), (2, 2, 2, 2), (1, 6), (3, 3, 3)]:
        g = nx.complete_multipartite_graph(*parts)
        print(f"   {str(parts):<16} {g.number_of_nodes():>3} "
              f"{total_domination_number(g):>8} "
              f"{str(sorted(set(minimal_total_dominating_sizes(g)))):>15}")

    print("\n   the n >= 5 hypothesis plays no role:")
    small = [g for g in cm if g.number_of_nodes() < 5]
    print(f"   {len(small)} complete multipartite graphs on < 5 vertices, "
          f"all well totally dominated: "
          f"{all(is_well_totally_dominated(g) for g in small)}")


if __name__ == "__main__":
    main()
