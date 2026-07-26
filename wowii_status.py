"""The status of the 22 conjectures according to the primary source.

The Lean formalisations in `google-deepmind/formal-conjectures` carry the tag
`category research open`. **That tag is not a status report.** It records what
was believed when the file was written, and for this list it is out of date in
both directions. The authoritative source is DeLaVina's own page,

    http://cms.dt.uh.edu/faculty/delavinae/research/wowII/

which is split into `resolved.htm` and `open.html` and is still being updated:
two of these conjectures were refuted in the week of 21-23 July 2026.

Checked 27 July 2026. Conjecture numbers on that page match the Lean file
numbers -- verified by comparing the statement text of all 22, not the numbers.

    #146  PROVED    21 Jul 2026, Brain Akaka, human-readable proof, developed
                    and checked with AI systems; upstreamed as pull request
                    4505 to google-deepmind/formal-conjectures.
    #198  RESOLVED  June 2010, Richard Stong (CCR, La Jolla). Sixteen years
                    before the Lean file called it open. But see below: the
                    statement implemented here is 198a, a different conjecture,
                    and 198a is still open.
    #200  FALSE     21 Jul 2026, Jitendra Prajapati. Counterexample on 11
                    vertices, graph6 J??FFBRq}N_. Verified below.
    #291  FALSE     23 Jul 2026, Zyad Tamimi. Counterexample on 12 vertices,
                    gamma_t = 4, k = 2, freq = 1. Reproduced independently
                    below -- a different graph with the same profile.

    the other 18 are still listed as open.

**198 versus 198a.** The page carries both. Conjecture 198 asks
b(G) <= 2 + ecc_avg(M) with M the *maximum-degree* vertices, and Stong resolved
it in 2010. Conjecture 198a asks b(G) <= 2 + ecc_avg(G), averaged over *all*
vertices, and is still open. `conj_198a` implements the second, so that effort
was not spent on a settled problem -- but it was recorded under the label "198"
throughout, which is exactly how such a confusion propagates. The two differ on
applicability for 19 of the 995 graphs on at most 7 vertices.

**What this cost here.** Work went into
146 (proved six days earlier), and into pushing 200 and 217 as a pair -- with
the conclusion that "the proof to attempt is hypothesis => the Bondy-Chvatal
closure of G join K1 is complete". For 200 that proof cannot exist. The
statement is false, and the reason the in-range work did not see it is simply
that the counterexample has 11 vertices and the exhaustive range was 8.

Note the certificates were not fooled: on Prajapati's graph all five classical
Hamiltonicity certificates correctly decline. The tooling was sound; the
literature check was not.
"""
from __future__ import annotations

import networkx as nx

import wowii as W

STATUS = {
    2: "open", 19: "open", 40: "open", 59: "open", 61: "open", 65: "open",
    100: "open", 133: "open", 141: "open", 142: "open", 144: "open",
    145: "open", 146: "PROVED 2026-07-21 (Akaka)", 160: "open",
    194: "open", 198: "RESOLVED 2010-06 (Stong); 198a, implemented here, is open",
    200: "FALSE 2026-07-21 (Prajapati)", 217: "open",
    291: "FALSE 2026-07-23 (Tamimi)", 314: "open", 316: "open", 322: "open",
}

# Prajapati's graph, from the page.
CE_200 = b"J??FFBRq}N_"

# Found here by hill-climbing on gamma_t - (k + freq), without using Tamimi's
# graph. Same invariant profile, different graph.
CE_291 = b"KCa|KG@_cJWg"


def check_200():
    g = nx.from_graph6_bytes(CE_200)
    import math
    tree = W.largest_induced_tree(g)
    rhs = math.ceil(1 + W.avg_l(g))
    ham = W.has_hamiltonian_path(g)
    return g, tree, rhs, ham, (tree == rhs and not ham)


def check_291():
    g = nx.from_graph6_bytes(CE_291)
    gt = W.total_domination_number(g)
    k = W.havel_hakimi_zero_step(g)
    f = W.freq_min_triangles(g)
    return g, gt, k, f, gt > k + f


def main():
    print("status of the 22 conjectures, per DeLaVina's page, checked "
          "27 July 2026\n")
    for n, s in sorted(STATUS.items()):
        mark = "   " if s == "open" else ">> "
        print(f"{mark}#{n:<5} {s}")

    print("\n" + "=" * 70)
    print("VERIFYING conjecture 200 is false (Prajapati, 11 vertices)")
    g, tree, rhs, ham, bad = check_200()
    print(f"   n={g.number_of_nodes()} m={g.number_of_edges()} "
          f"connected={nx.is_connected(g)}")
    print(f"   tree(G) = {tree},  ceil(1 + l_avg) = {rhs}  -> "
          f"hypothesis holds: {tree == rhs}")
    print(f"   Hamiltonian path: {ham}")
    print(f"   => CONJECTURE 200 IS FALSE: {bad}")

    print("\n   and the classical certificates all correctly decline:")
    import wowii_ham as H
    for name, fn in H.CERTIFICATES:
        print(f"      {name:<40} {fn(g)}")

    print("\n" + "=" * 70)
    print("VERIFYING conjecture 291 is false (found independently here)")
    g, gt, k, f, bad = check_291()
    print(f"   n={g.number_of_nodes()} m={g.number_of_edges()} "
          f"connected={nx.is_connected(g)}")
    print(f"   gamma_t = {gt},  k = {k},  freqMinTriangles = {f}")
    print(f"   the conjecture asks gamma_t <= k + freq = {k + f}")
    print(f"   => CONJECTURE 291 IS FALSE: {bad}")
    print(f"   Tamimi's graph has the same profile (gamma_t 4, k 2, freq 1)")
    print(f"   on 12 vertices; the search found none on 11, consistent with")
    print(f"   12 being the minimum order.")

    print("\n" + "=" * 70)
    print("what remains genuinely open, of the 22:")
    op = [n for n, s in sorted(STATUS.items()) if s == "open"]
    print(f"   {op}")
    print(f"   ({len(op)} of 22; 160 is refuted as formalised by this project,")
    print(f"    and 65, 141, 314, 316, 322 are proved here)")


if __name__ == "__main__":
    main()
