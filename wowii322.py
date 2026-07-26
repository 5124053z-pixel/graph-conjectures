"""Written on the Wall II, Conjecture 322 — resolved, trivially.

    WOWII 322. Let G be a simple connected graph on n >= 5 vertices. If
    max_v l(v) <= 1, where l(v) = alpha(G[N(v)]) is the independence number of
    the neighbourhood of v, then G is well totally dominated.

Listed as open, and formalised as open in Google DeepMind's `formal-conjectures`
repository (`WrittenOnTheWallII/GraphConjecture322.lean`, tagged
`category research open`).

    THEOREM. The conjecture is true, and the hypothesis n >= 5 is unnecessary:
    it holds for every connected graph on at least two vertices.

Proof, in two steps.

**Step 1: the hypothesis forces G to be complete.**
`l(v) <= 1` says the induced subgraph on N(v) has no two non-adjacent vertices,
i.e. N(v) is a clique. Suppose u ~ v ~ w with u != w. Then u, w in N(v), so
u ~ w. Adjacency is therefore transitive, so "u = v or u ~ v" is an equivalence
relation and every connected component is a clique. G is connected, hence
G = K_n.

**Step 2: K_n is well totally dominated.**
In K_n with n >= 2 every 2-element subset {u, v} is a total dominating set: any
vertex w has a neighbour in it (if w != u then u is a neighbour of w; if w = u
then v is). A singleton {u} is not, since u itself has no neighbour in {u}. And
any set of size >= 3 properly contains a total dominating pair, so it is not
minimal. The minimal total dominating sets are therefore exactly the pairs, all
of size 2, and G is well totally dominated. QED

**On what this is worth.** Very little as mathematics -- it is two lines, and
the reason nobody had written them is presumably that nobody looked. What it
illustrates is a property of automatically generated conjectures: the hypothesis
`max_v l(v) <= 1` *reads* like a mild local condition and is in fact the
strongest possible one, collapsing the statement to a single family. Among all
995 connected graphs on at most 7 vertices, exactly 6 satisfy it, and they are
K_2 through K_7.

That is the kind of thing worth checking before spending a search budget on a
conjecture.
"""
from __future__ import annotations

import networkx as nx

from wowii import (local_independence, is_well_totally_dominated,
                   minimal_total_dominating_sizes)


def hypothesis_holds(g):
    return all(local_independence(g, v) <= 1 for v in g)


def main():
    print("WOWII 322: connected, n >= 5, max_v l(v) <= 1  =>  "
          "well totally dominated\n")

    print("step 1 -- the hypothesis forces G to be complete.")
    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= 7 and nx.is_connected(g)]
    sat = [g for g in graphs if hypothesis_holds(g)]
    n = len(graphs)
    complete = all(2 * g.number_of_edges()
                   == g.number_of_nodes() * (g.number_of_nodes() - 1)
                   for g in sat)
    print(f"   of the {n} connected graphs on 2..7 vertices, {len(sat)} satisfy it")
    print(f"   orders: {sorted(g.number_of_nodes() for g in sat)}")
    print(f"   every one of them complete: {complete}")

    print("\nstep 2 -- K_n is well totally dominated, for every n >= 2.")
    print(f"   {'graph':<8} {'minimal total dominating set sizes':<36} WTD")
    ok = True
    for k in range(2, 10):
        K = nx.complete_graph(k)
        sizes = sorted(minimal_total_dominating_sizes(K))
        wtd = is_well_totally_dominated(K)
        ok &= wtd and sizes == [2]
        print(f"   K_{k:<6} {str(sizes):<36} {wtd}")

    print(f"\nboth steps check out: {complete and ok}")
    print("the conjecture holds, and n >= 5 was never needed")


if __name__ == "__main__":
    main()
