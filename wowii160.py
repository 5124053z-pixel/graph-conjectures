"""Written on the Wall II, Conjecture 160 — the Lean file is wrong, the
conjecture is not refuted.

    THE HEADLINE, stated first because an earlier version of this file buried
    it. Conjecture 160 is NOT refuted by anything here. What is refuted is the
    transcription in `GraphConjecture160.lean`, and the evidence is now
    decisive rather than a single counterexample:

        c_C4 read as a COUNT of induced 4-cycles (the Lean file):
            8,985 violations among the 12,112 connected graphs on <= 8 vertices
        c_C4 read as an INDICATOR, 0 if G has a C4 and 1 otherwise
            -- exactly the convention conjecture 133 uses for the same symbol:
            0 violations

    A statement that fails on three quarters of all small graphs is not a
    conjecture anyone would publish, and the same symbol is defined as an
    indicator ten conjectures earlier on the same list. So the count reading is
    a transcription error, and under the indicator reading conjecture 160 holds
    everywhere it has been checked and remains open.

`formal-conjectures/FormalConjectures/WrittenOnTheWallII/GraphConjecture160.lean`
states, tagged `@[category research open]`:

    L_s(G)  >=  max_v l(v)  +  max_v T(v) * c_{C4}(G)

for every simple connected G, where L_s is the maximum number of leaves over all
spanning trees, l(v) is the independence number of N(v), T(v) is the number of
triangles containing v, and c_{C4}(G) is the **number** of induced 4-cycles
(`countInducedC4`, which counts unordered 4-sets inducing a C4).

    COUNTEREXAMPLE. G = K_5 minus a perfect matching on four of its vertices:
    V = {0,1,2,3,4}, E = {01, 03, 04, 12, 14, 23, 24, 34}
    (equivalently, K_5 with the edges 02 and 13 removed).

Every quantity is checkable by hand.

* Degrees are 3,3,3,3,4 — vertex 4 is adjacent to all others.
* `N(4) = {0,1,2,3}` induces exactly the 4-cycle 0–1–2–3–0, since 02 and 13 are
  the removed edges. So `l(4) = α(C_4) = 2`. Every other neighbourhood induces a
  path on three vertices, also with independence number 2. **max_v l(v) = 2.**
* The triangles are {4,0,1}, {4,1,2}, {4,2,3}, {4,3,0}; none avoids vertex 4,
  because {0,1,2,3} induces a C_4 and has no triangle. **max_v T(v) = 4.**
* The only 4-set inducing a C_4 is {0,1,2,3}: every other 4-set contains vertex
  4, which is adjacent to the other three, so it induces at least 5 edges.
  **c_{C4}(G) = 1.**
* A spanning tree on 5 vertices has 4 edges and so at most 4 leaves; the star
  centred at vertex 4 attains it. **L_s(G) = 4.**

Therefore the right-hand side is `2 + 4 * 1 = 6 > 4 = L_s(G)`. QED

**What this does and does not show.** It refutes the statement as formalised. It
does not settle DeLaVina's conjecture 160, because the formalisation may not
render it faithfully — and there is direct evidence of that risk in the same
directory:

* conjecture 133 defines `cC4` as an **indicator**, `if hasC4 then 0 else 1`,
  while conjecture 160 uses the same name for a **count**;
* conjecture 100's doc comment says `diam(Gᶜ)` where its Lean statement uses
  `degreeL2Norm Gᶜ`.

So the likely reading is that `c_{C4}` in the original is an indicator too — and
it does rescue the statement. **An earlier version of this file said otherwise,
and was wrong because of a mislabelled line.** It tested "0 if C4, conjecture
133's convention" while actually computing the indicator of an *induced*
4-cycle; conjecture 133 uses `hasC4`, a 4-cycle as a **subgraph**, and those two
tests differ on 310 of the 995 connected graphs on at most 7 vertices. Under the
real convention the statement has **no violation anywhere it has been checked**.

The useful output here is a bug report against the formalisation, not a
mathematical claim.
"""
from __future__ import annotations

import networkx as nx

from wowii import (has_c4_subgraph, local_independence, max_triangles_at_vertex,
                   count_induced_c4, max_leaves_spanning_tree)

COUNTEREXAMPLE = nx.Graph([(0, 1), (0, 3), (0, 4), (1, 2),
                           (1, 4), (2, 3), (2, 4), (3, 4)])


def readings(g):
    """The four ways c_{C4} can be read.

    The third line is the one that matters and it was **mislabelled** in an
    earlier version of this file: it said "conjecture 133's convention" while
    computing the indicator of an *induced* 4-cycle. Conjecture 133 uses
    `hasC4`, a 4-cycle as a **subgraph**, which is a different test on 310 of
    the 995 connected graphs on at most 7 vertices. Under the real convention
    the statement never fails; under the mislabelled one it fails often, and
    that is how this file came to claim that no reading rescues it."""
    n4 = count_induced_c4(g)
    sub = has_c4_subgraph(g)
    return {
        "count of induced C4 (as formalised)": n4,
        "indicator, 0 if INDUCED C4 else 1": 0 if n4 else 1,
        "indicator, 1 if induced C4 else 0": 1 if n4 else 0,
        "indicator, 0 if C4 SUBGRAPH else 1  (conjecture 133's actual convention)":
            0 if sub else 1,
    }


def main():
    g = COUNTEREXAMPLE
    print("WOWII 160:  Ls(G) >= max_v l(v) + max_v T(v) * c_C4(G)\n")
    print(f"counterexample: K_5 minus the edges 02 and 13")
    print(f"   V = {sorted(g)}")
    print(f"   E = {sorted(g.edges())}")
    print(f"   degrees            {dict(sorted(g.degree()))}")
    ml = max(local_independence(g, v) for v in g)
    mt = max_triangles_at_vertex(g)
    c4 = count_induced_c4(g)
    ls = max_leaves_spanning_tree(g)
    print(f"   max_v l(v)         {ml}")
    print(f"   max_v T(v)         {mt}")
    print(f"   c_C4(G)            {c4}")
    print(f"   Ls(G)              {ls}")
    print(f"\n   right-hand side    {ml} + {mt} * {c4} = {ml + mt * c4}")
    print(f"   left-hand side     {ls}")
    print(f"   {ml + mt * c4} > {ls}, so the statement fails."
          f"  {'CONFIRMED' if ml + mt * c4 > ls else 'DOES NOT FAIL'}")

    print("\nsmallest counterexample per reading,")
    print("over all connected graphs on 2..6 vertices:")
    graphs = [x for x in nx.graph_atlas_g()
              if 2 <= x.number_of_nodes() <= 6 and nx.is_connected(x)]
    for name in readings(g):
        bad = []
        for x in graphs:
            ml = max(local_independence(x, v) for v in x)
            val = readings(x)[name]
            if ml + max_triangles_at_vertex(x) * val > max_leaves_spanning_tree(x):
                bad.append(x.number_of_nodes())
        print(f"   {name:<50} {len(bad):>3} violations, smallest n = "
              f"{min(bad) if bad else '-'}")



def compare_readings(nmax=8):
    """The decisive test: which reading of c_C4 makes the statement plausible."""
    import networkx as nx
    from wowii import (largest_induced_bipartite, local_independence,
                       max_triangles_at_vertex, count_induced_c4,
                       has_c4_subgraph, max_leaves_via_cds)
    if nmax <= 7:
        gs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= nmax and nx.is_connected(g)]
    else:
        from allgraphs import connected_graphs
        gs = connected_graphs(nmax, verbose=False)
    ml = lambda g: max(local_independence(g, v) for v in g)
    readings = [
        ("count of induced C4  (the Lean file)", count_induced_c4),
        ("indicator, 0 if G has a C4 else 1  (as in conjecture 133)",
         lambda g: 0 if has_c4_subgraph(g) else 1),
        ("indicator on induced C4",
         lambda g: 0 if count_induced_c4(g) > 0 else 1),
    ]
    print()
    print(f"the two readings of c_C4, over the {len(gs)} connected graphs "
          f"on 2..{nmax}")
    for name, c in readings:
        bad = sum(1 for g in gs
                  if max_leaves_via_cds(g) < ml(g) + max_triangles_at_vertex(g) * c(g))
        verdict = "HOLDS EVERYWHERE" if not bad else f"{bad} violations"
        print(f"   {name:<58} {verdict}")
    print()
    print("   A statement failing on three quarters of all small graphs is not")
    print("   one anybody lists as open for twenty years. The count reading is a")
    print("   transcription error; conjecture 160 itself is NOT refuted here.")


if __name__ == "__main__":
    main()
    compare_readings(8)
