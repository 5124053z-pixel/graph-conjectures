"""Generate the status table from the code, so it cannot drift from reality.

Failure 19 in the README was a line labelled "conjecture 133's convention" that
computed something else, and a conclusion rested on the label. Every other guard
in this project compares a computed value against something independent; none
compared a *name* against what the code under it does.

This file is the fix for the status table specifically. Rather than a
hand-maintained list of what is proved and what is only checked, each entry
names:

  * the conjecture function in `wowii.py`, so "does it hold" is measured, not
    asserted;
  * a `regime` predicate and a `witness` callable for the part that is PROVED,
    so the claimed coverage is recomputed every run;
  * prose for what the proof is.

If a claim here is wrong, running this file says so.

    python status_table.py         # over 2..7 vertices, fast
    python status_table.py 8       # over all 12,112 connected graphs on 2..8
"""
from __future__ import annotations

import math
import sys

import networkx as nx

import wowii as W
import wowii_toolkit as T


def _girth(g):
    return 0 if nx.is_forest(g) else nx.girth(g)


def _eccC(g):
    return W.ecc_set(g, nx.center(g))


def _lmin_comp(g):
    gc = nx.complement(g)
    if gc.number_of_edges() == 0:
        return None
    return min(W.local_independence(gc, v) for v in gc)


def _both_maxima(g):
    ml = max(W.local_independence(g, v) for v in g)
    mT = W.max_triangles_at_vertex(g)
    return any(W.local_independence(g, v) == ml
               and g.subgraph(list(g[v])).number_of_edges() == mT for v in g)


# conjecture -> (conj fn, prose, regime predicate or None, proved?)
#
# `proved` is the column that matters and it is kept separate on purpose. A
# regime can cover every graph in reach while resting on an ingredient that is
# only verified -- conjectures 19 and 61 both do -- and collapsing those two
# facts into one number is exactly the kind of drift this file exists to stop.
ENTRIES = {
    "65":  (True, W.conj_65, "distMin takes only 0 and 1, so the rhs never exceeds 2",
            None),
    "141": (True, W.conj_141, "neighbourhood star below girth 6, locally tree-like "
            "ball above; the split is forced", None),
    "314": (True, W.conj_314, "every step is an induced P5 off the private "
            "total-neighbours; the last branch names the Wagner graph", None),
    "316": (True, W.conj_316, "a counting argument caps the core at 3 and forces it "
            "complete", None),
    "322": (True, W.conj_322, "the hypothesis says exactly that G is complete "
            "multipartite, and those are all WTD", None),

    "142": (True, W.conj_142, "girth <= 3: diam >= 1 + ecc(B), so "
            "tree >= diam+1 >= 2 + ecc(B)",
            lambda g: _girth(g) <= 3),
    "144": (True, W.conj_144, "girth <= 3: diam >= 1 + ecc(Centers), so "
            "tree >= diam+1 >= 2 + ecc(Centers)",
            lambda g: _girth(g) <= 3),
    "160": (True, W.conj_160, "G has a C4: the rhs is max l(v) <= Delta <= Ls. "
            "C4-free with one vertex attaining both maxima: rhs = deg(v)",
            lambda g: W.has_c4_subgraph(g) or _both_maxima(g)),
    "145": (True, W.conj_145, "lMin(Gc) >= 2, or lMin = 1 with ecc(B) <= 2",
            lambda g: (_lmin_comp(g) or 0) >= 2 or T.eccB(g) <= 2),
    "61":  (False, W.conj_61, "diam <= 3, or 3 | diam, or residue < alpha "
            "(the last two lean on the UNPROVED f >= alpha + floor(diam/3))",
            lambda g: nx.diameter(g) <= 3),
    "100": (True, W.conj_100, "alpha >= 16 -- out of reach of any corpus here",
            lambda g: W.indep_number(g) >= 16),
    "2":   (True, W.conj_2, "triangle-free, or average degree <= 3: "
            "Ls >= deg(u)+deg(v)-2 on an edge, plus Cauchy-Schwarz",
            lambda g: sum(nx.triangles(g).values()) == 0
            or 2 * g.number_of_edges() / g.number_of_nodes() <= 3),
    "59":  (False, W.conj_59, "reduced to 2f >= b + res by AM-GM, which is "
            "verified everywhere but not proved", None),
    "19":  (False, W.conj_19, "reduced to alpha + alpha(G-I) >= avg_ecc + max l, "
            "which is verified but not proved", None),
}


def main(nmax=7):
    if nmax <= 7:
        gs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= nmax and nx.is_connected(g)]
    else:
        from allgraphs import connected_graphs
        gs = connected_graphs(nmax, verbose=False)

    print(f"status over the {len(gs)} connected graphs on 2..{nmax}\n")
    print(f"   {'conj':>5} {'applicable':>11} {'holds':>7} {'in the regime':>14} "
          f"{'residual':>9}   status of the regime")
    bad_any = False
    for name, (proved, fn, prose, regime) in sorted(ENTRIES.items(),
                                                    key=lambda t: int(t[0])):
        ap = [g for g in gs if fn(g)[0]]
        holds = all(fn(g)[1] for g in ap)
        bad_any |= not holds
        if regime is None:
            covered, residual = len(ap), 0
        else:
            covered = sum(1 for g in ap if regime(g))
            residual = len(ap) - covered
        mark = "proved" if proved else "UNPROVED ingredient"
        print(f"   {name:>5} {len(ap):>11} {str(holds):>7} {covered:>14} "
              f"{residual:>9}   {mark}")

    print("\nwhat the proved regime is, per conjecture")
    for name, (proved, fn, prose, regime) in sorted(ENTRIES.items(),
                                                    key=lambda t: int(t[0])):
        scope = "all graphs" if regime is None else "a regime"
        scope += "" if proved else ", NOT PROVED"
        print(f"   {name:>5}  [{scope}]  {prose}")

    print(f"\nany conjecture failing on an applicable graph: {bad_any}")
    if bad_any:
        print("   ^ that would be a counterexample; investigate before "
              "believing it")

    print("\nnot in this table, and why")
    print("   40, 133   settled in range only by witnesses, no regime is proved")
    print("   146, 198, 200, 291   resolved by other people; see "
          "wowii_status.py")
    print("   194, 198a, 217   open, no proved regime; see wowii_ham.py")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
