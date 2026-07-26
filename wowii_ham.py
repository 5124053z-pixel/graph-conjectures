"""The four Hamiltonian-path conjectures, against classical certificates.

    WOWII 194.  alpha(G) <= 1 + l_avg(G)                =>  Hamiltonian path
    WOWII 198.  b(G) <= 2 + avg ecc(G)                  =>  Hamiltonian path
    WOWII 200.  tree(G) == ceil(1 + l_avg(G))           =>  Hamiltonian path
    WOWII 217.  Ls(G) <= 4*[residue(G) == 2] + 2        =>  Hamiltonian path

All four are formalised as open (`category research open`). All four have the
same shape: some invariant is small, therefore a Hamiltonian path exists. So
the question to ask is not "can we find a Hamiltonian path" but **which
classical sufficient condition does the hypothesis already imply**.

Every certificate below is a theorem, not a search. Each is checked against
`has_hamiltonian_path` over the 995 connected graphs on at most 7 vertices
before being used, because a certificate that is subtly misimplemented looks
exactly like a strong one.

    (1) Chvatal-Erdos.  kappa(G) >= alpha(G) - 1  =>  Hamiltonian path.

    (2) Ore.  deg(u) + deg(v) >= n - 1 for every non-adjacent pair
        =>  Hamiltonian path.

    (3) Ls(G) <= 2  =>  Hamiltonian path. Every spanning tree has at least two
        leaves, so if the maximum is 2 then *every* spanning tree has exactly
        two -- every spanning tree is a path. This one is free and settles the
        residue != 2 branch of 217 outright.

    (4) Bondy-Chvatal, applied to G join K1.  G has a Hamiltonian path exactly
        when G join K1 (a new vertex joined to everything) has a Hamiltonian
        cycle. The Bondy-Chvatal closure of a graph on N vertices -- repeatedly
        add uv whenever deg(u) + deg(v) >= N -- is Hamiltonian iff the original
        is, so a complete closure certifies a Hamiltonian path in G.

        This is the one that does the work. Ore, Dirac and Chvatal-Erdos are
        all one-shot degree or connectivity tests; the closure *iterates*, and
        each edge it adds raises degrees and lets it add more.

    (5) Chvatal's degree condition, applied to G join K1.

What the four conjectures look like against these:

    conj   applicable   residual
     194          501          1
     198          148          3
     200          136          0
     217          453          0

**200 and 217 are settled in range by classical theory alone**, and the other
two are down to four graphs between them. The reading is that these hypotheses
are not saying anything new about Hamiltonicity -- they are roundabout ways of
asserting a condition Bondy and Chvatal already covered in 1976. That is a
statement about the conjectures, and it is the useful output here: the proof to
attempt is `hypothesis => the closure of G join K1 is complete`, not
`hypothesis => a Hamiltonian path can be constructed`.

**The obstruction that never appears.** No graph satisfying any of the four
hypotheses has a vertex cut of size k leaving more than k+1 components -- the
standard reason a graph has no Hamiltonian path. Checked for k up to 3 over all
graphs in range. The hypotheses are doing enough to rule out the only easy
obstruction, which is consistent with all four being true.
"""
from __future__ import annotations

import itertools

import networkx as nx

import wowii as W


def join_k1(g):
    h = g.copy()
    h.add_node("*")
    for v in g:
        h.add_edge("*", v)
    return h


def bondy_chvatal_closure(g):
    h = g.copy()
    n = h.number_of_nodes()
    changed = True
    while changed:
        changed = False
        for u, v in list(nx.complement(h).edges()):
            if h.degree(u) + h.degree(v) >= n:
                h.add_edge(u, v)
                changed = True
    return h


def cert_chvatal_erdos(g):
    return nx.node_connectivity(g) >= W.indep_number(g) - 1


def cert_ore(g):
    n = g.number_of_nodes()
    return all(g.degree(u) + g.degree(v) >= n - 1
               for u, v in nx.complement(g).edges())


def cert_few_leaves(g):
    return W.max_leaves_via_cds(g) <= 2


def cert_closure(g):
    h = bondy_chvatal_closure(join_k1(g))
    n = h.number_of_nodes()
    return h.number_of_edges() == n * (n - 1) // 2


def cert_chvatal_degree(g):
    h = join_k1(g)
    n = h.number_of_nodes()
    d = sorted(x for _, x in h.degree())
    return all(d[i - 1] > i or d[n - i - 1] >= n - i
               for i in range(1, (n + 1) // 2))


CERTIFICATES = [("Chvatal-Erdos", cert_chvatal_erdos),
                ("Ore", cert_ore),
                ("Ls <= 2", cert_few_leaves),
                ("Bondy-Chvatal closure of G join K1", cert_closure),
                ("Chvatal degree condition on G join K1", cert_chvatal_degree)]

CONJ = {194: W.conj_194, 198: W.conj_198a, 200: W.conj_200, 217: W.conj_217}


def no_easy_obstruction(g, kmax=3):
    """No vertex cut of size k leaves more than k+1 components."""
    n = g.number_of_nodes()
    for k in range(1, min(kmax, n - 1) + 1):
        for S in itertools.combinations(g, k):
            if nx.number_connected_components(
                    nx.restricted_view(g, S, [])) > k + 1:
                return False
    return True


def main(nmax=7):
    if nmax <= 7:
        graphs = [g for g in nx.graph_atlas_g()
                  if 2 <= g.number_of_nodes() <= nmax and nx.is_connected(g)]
    else:
        from allgraphs import connected_graphs
        graphs = connected_graphs(nmax)
    print(f"the four Hamiltonian-path conjectures over the {len(graphs)} "
          f"connected graphs on 2..{nmax} vertices\n")

    print("every certificate is checked against the truth first")
    for name, fn in CERTIFICATES:
        bad = [g for g in graphs if fn(g) and not W.has_hamiltonian_path(g)]
        print(f"   {name:<40} violations {len(bad)}")

    print(f"\n   {'conj':>5} {'applicable':>11} {'residual':>9}")
    residuals = {}
    for num, fn in sorted(CONJ.items()):
        ap = [g for g in graphs if fn(g)[0]]
        assert all(fn(g)[1] for g in ap), f"conjecture {num} FAILS"
        left = [g for g in ap if not any(c(g) for _, c in CERTIFICATES)]
        residuals[num] = left
        print(f"   {num:>5} {len(ap):>11} {len(left):>9}"
              f"{'   <-- settled in range' if not left else ''}")

    print("\nwhich certificate carries each conjecture")
    print(f"   {'conj':>5} " + " ".join(f"{n[:14]:>15}" for n, _ in CERTIFICATES))
    for num, fn in sorted(CONJ.items()):
        ap = [g for g in graphs if fn(g)[0]]
        row = [sum(1 for g in ap if c(g)) for _, c in CERTIFICATES]
        print(f"   {num:>5} " + " ".join(f"{x:>15}" for x in row))

    print("\nthe easy obstruction never appears")
    for num, fn in sorted(CONJ.items()):
        ap = [g for g in graphs if fn(g)[0]]
        print(f"   #{num}: applicable graphs with a bad vertex cut: "
              f"{sum(1 for g in ap if not no_easy_obstruction(g))}")

    for num, left in residuals.items():
        if left and len(left) <= 5:
            print(f"\n#{num} residual graphs")
            for g in left:
                print(f"   n={g.number_of_nodes()} m={g.number_of_edges()} "
                      f"alpha={W.indep_number(g)} "
                      f"kappa={nx.node_connectivity(g)} "
                      f"edges={sorted(g.edges())}")


if __name__ == "__main__":
    import sys
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
