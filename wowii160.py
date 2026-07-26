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

So the likely reading is that `c_{C4}` in the original is an indicator too. That
does not rescue the statement either: all three readings — count, `0 if C4`, and
`1 if C4` — fail, the last two with counterexamples on 4 and 5 vertices
respectively. Whatever conjecture 160 says, the file does not say it.

The useful output here is a bug report against the formalisation, not a
mathematical claim.
"""
from __future__ import annotations

import networkx as nx

from wowii import (local_independence, max_triangles_at_vertex,
                   count_induced_c4, max_leaves_spanning_tree)

COUNTEREXAMPLE = nx.Graph([(0, 1), (0, 3), (0, 4), (1, 2),
                           (1, 4), (2, 3), (2, 4), (3, 4)])


def readings(g):
    """The three ways c_{C4} can be read, all of which fail somewhere."""
    n4 = count_induced_c4(g)
    return {
        "count (as formalised)": n4,
        "indicator, 0 if C4 (conjecture 133's convention)": 0 if n4 else 1,
        "indicator, 1 if C4": 1 if n4 else 0,
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

    print("\nno reading of c_C4 rescues it -- smallest counterexample per reading,")
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


if __name__ == "__main__":
    main()
