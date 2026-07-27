r"""WOWII 59:  f(G) >= ceil( sqrt( res(G) * b(G) ) ).

    f    = forest number, the order of a largest induced forest
    res  = Havel-Hakimi residue
    b    = bipartite number, the order of a largest induced bipartite subgraph

Everything below is labelled PROVED or CHECKED and the labels are not
decorative. PROVED means there is an argument in the docstring that works for
every graph. CHECKED means a computation over a finite corpus and nothing more;
conjecture 200 of this same list held over all 12,112 connected graphs on at
most 8 vertices and is false at 11, so CHECKED is worth exactly the size of the
corpus.

-------------------------------------------------------------------------------
PROVED
-------------------------------------------------------------------------------

P1.  f(G) <= b(G).
     A forest is bipartite, so a largest induced forest is an induced bipartite
     subgraph.

P2.  f(G) >= alpha(G) + 1  for connected G on n >= 2 vertices.
     Let I be a maximum independent set. Since G is connected and n >= 2 there
     is a vertex v outside I; G[I u {v}] is a star centred at v, hence a forest.
     (Already used in this project.)

P3.  b(G) <= 2*alpha(G).
     The two colour classes of a largest induced bipartite subgraph are
     independent sets of G, each of size at most alpha.

P4.  f(G) >= ceil(b(G)/2) + 1.
     P3 gives alpha >= b/2, and alpha is an integer, so alpha >= ceil(b/2);
     apply P2. This is DeLaVina's note of 6 March 2004 under conjecture 40; the
     two-line derivation from P2 and P3 is the point, not the statement.

P5.  res(G) <= alpha(G).
     Favaron, Maheo and Sacle, *On the residue of a graph*, JGT 15 (1991).
     Cited, not reproved here.

P6.  res(G) <= n - ceil(m/Delta).
     Each Havel-Hakimi step deletes one entry of current value d and removes
     exactly d from the remaining degree sum's half, so the degrees removed over
     all k steps sum to m; every step removes at most Delta, so k >= m/Delta,
     and res = n - k.
     (Valid but weak -- it never decides anything the other bounds do not.)

L2.  If some maximum independent set I admits two distinct vertices a, b outside
     I with G[I u {a,b}] acyclic, then f(G) >= alpha(G) + 2.
     Immediate: the set is a witness. What makes it usable is that the condition
     is a search over O(n^2) pairs per maximum independent set, not over subsets.
     Explicitly, G[I u {a,b}] is acyclic iff  |N(a) n N(b) n I| <= 1  and, when
     ab is an edge, |N(a) n N(b) n I| = 0 -- one common neighbour plus the edge
     ab is a triangle, two common neighbours are a 4-cycle, and no other cycle
     can exist because I is independent.

T1 (case i).   res*b <= (alpha+1)^2  =>  conjecture 59.
     By P2, f >= alpha+1, so f^2 >= (alpha+1)^2 >= res*b; take square roots and
     use that f is an integer to absorb the ceiling.

T1'.  res <= 2  =>  conjecture 59.
     res*b <= 2*b <= 4*alpha <= (alpha+1)^2, since (alpha-1)^2 >= 0. So T1
     applies. This is the whole of the "free" range and it is not empty:
     4,382 of the 12,112 graphs in range have res <= 2.

T2 (case ii).  b <= f+1  =>  conjecture 59.
     res <= alpha <= f-1 by P5 and P2, so res*b <= (f-1)(f+1) = f^2 - 1 < f^2.

T3 (case iii). L2 applies and res*b <= (alpha+2)^2  =>  conjecture 59.

T4 (the AM-GM reduction).  2*f >= b + res  =>  conjecture 59.
     sqrt(res*b) <= (res+b)/2 <= f by AM-GM, and f is an integer, so
     f >= ceil(sqrt(res*b)). The converse fails, so this is a strictly stronger
     statement; see CHECKED C2.

T5.  res <= 2  =>  2*f >= b + res.
     2f >= 2*alpha + 2 >= b + 2 >= b + res, by P2 and P3.

T6.  f >= alpha + res/2  =>  2*f >= b + res.
     b <= 2*alpha <= 2f - res by P3. The hypothesis is FALSE in general -- see
     FAILED below -- so this is a dead end, recorded because it is the obvious
     first thing to try.

-------------------------------------------------------------------------------
CHECKED (finite corpus, no proof)
-------------------------------------------------------------------------------

C1.  Conjecture 59 holds on all 12,112 connected graphs on 2..8 vertices.
     Already in the repository; re-verified here with an independent bitmask
     implementation that is first validated against wowii.py on all 995 graphs
     on at most 7 vertices.

C2.  **2*f >= b + res holds on all 12,112 of them, with zero violations.**
     By T4 this is strictly stronger than conjecture 59, and it is sharp: 196
     of the 12,112 attain equality, and equality holds for K_{a,a} for every a
     (f = a+1, b = 2a, res = 2), so it cannot be improved to 2f >= b + res + 1.
     Replacing res by alpha -- i.e. 2f >= b + alpha -- is FALSE (282
     violations), so the strengthening genuinely needs the residue.

C3.  Cases (i) and (iii) alone cover all 12,112 graphs: residual 0, for
     conjecture 59 *and* for the stronger C2. Two cases, no girth or diameter
     split. (Case (ii) is then redundant in range but is kept because it is the
     cheapest and it is what makes the b ~ f regime obviously free.)

C4.  **The case analysis does not scale, and the point where it stops is
     computable in advance.** Case (iii) needs res*b <= (alpha+2)^2 while P3
     and P5 allow res*b to be as large as 2*alpha^2; 2*alpha^2 <= (alpha+2)^2
     exactly when alpha <= 4. Every one of the 14 residual graphs on <= 8
     vertices has alpha in {3,4} -- at the very edge -- and a hill-climb finds
     an alpha = 5 graph on 11 vertices where the structural bounds fail:

         JDAp@_PEF??   n=11, m=15, alpha = res = 5, b = 10
         best structural lower bound 7, ceil(sqrt(res*b)) = 8, true f = 9.

     So the conjecture still holds there; the *bounds* do not. Same signature as
     conjectures 145 and 146, whose residuals were 0 at n <= 7 and 75 and 96 at
     n <= 8. This is the honest verdict on the whole approach.

C5.  Counterexample hunts above the exhaustive range found nothing:
       - wowii_hunt.py on 59, n = 9..13: best margin 0 at every order.
       - the AM-GM strengthening 2f >= b + res, n = 9..14: best margin 0.
       - a hunt aimed at the crux regime res = alpha, b = 2*alpha, n = 9..13:
         reaches that regime at every order, and there the slack *grows*
         (1, 1, 2, 2, 3 for alpha = 3, 3, 4, 4, 5).
       - 2,989 random connected graphs on 9..12 vertices, edge probability
         uniform on [0.12, 0.7]: no violation of either statement, and both
         reach slack 0 at every order.
     Margin 0 means equality was reached, so the search was in the right
     neighbourhood. It is still weak evidence: the same method failed to
     rediscover conjecture 200's counterexample in ten minutes.

C6.  Where the two statements are tight is *different*, and that is the most
     useful thing found here. Conjecture 59 has slack growing without bound on
     K_{a,a} (f = a+1 against ceil(2*sqrt(a))), while the AM-GM strengthening is
     exactly tight there. And K_{a,a} has res = 2, which is inside T1', the
     regime that is free. **The strengthening is tight where it is easy and slack
     where it is hard** -- in the crux regime res = alpha, b = 2*alpha, the
     observed slack grows with alpha rather than shrinking.

-------------------------------------------------------------------------------
FAILED (recorded so they are not tried again)
-------------------------------------------------------------------------------

F1.  f >= alpha + res/2.   FALSE, 1,093 violations. Smallest: the star K_{1,3}
     (f = 4, alpha = 3, res = 3). Would have given the AM-GM form outright by T6.

F2.  f >= (alpha + b)/2.   FALSE, 282 violations, smallest on 6 vertices
     (graph6 EEz_: f = 4, alpha = 3, b = 6). So res cannot be replaced by alpha.

F3.  res*b <= (alpha+1)^2.  FALSE, 1,899 violations. P2 alone is not enough;
     something has to prove f >= alpha + 2.

F4.  b <= 2*res + 2.  FALSE, 159 violations.

F5.  f >= b - res + 1.  FALSE, 197 violations.

F6.  "f >= alpha + f(G[U_I]) where U_I is the set of vertices outside I with at
     most one neighbour in I."  FALSE, 957 violations. The tempting
     generalisation of the caterpillar bound: it is true when the added set is
     independent (two vertices with one neighbour each in I need two shared
     neighbours to close a 4-cycle) but false as soon as the added set has an
     edge -- adjacent u, v sharing their single I-neighbour x close the triangle
     uvx. The corrected version needs the neighbour map to be injective, which
     is a matching condition and no longer a clean bound.

F7.  Every upper bound on res tried other than P5 and P6 is either false
     (res <= n - Delta + 1, res + b <= n + 2, res <= 2n - m - 1) or useless.
     Bounding res from above is the crux and nothing here does it: residue
     depends only on the degree sequence, while b and f do not, so any proof has
     to convert a degree-sequence fact into a structural one.

-------------------------------------------------------------------------------
WHAT WOULD ACTUALLY SETTLE IT
-------------------------------------------------------------------------------

By T4, P3 and P5 the whole conjecture reduces to the regime where res is close
to alpha *and* b is close to 2*alpha, and there the sharp question is:

    Does  res(G) = alpha(G)  and  b(G) = 2*alpha(G)  force  f(G) >= ceil(sqrt(2)*alpha)?

That is the only case the free bounds cannot reach, it is where the 14 residual
graphs live, and the crux hunt in C5 says the true slack there grows rather than
closing. Nothing in this file proves it.
"""
from __future__ import annotations

import itertools
import math
import sys

import networkx as nx

import wowii as W


# ---------------------------------------------------------------------------
# fast exact invariants (validated against wowii.py before use)
# ---------------------------------------------------------------------------

def _masks(g):
    nodes = sorted(g)
    idx = {v: i for i, v in enumerate(nodes)}
    adj = [0] * len(nodes)
    for u, v in g.edges():
        adj[idx[u]] |= 1 << idx[v]
        adj[idx[v]] |= 1 << idx[u]
    return len(nodes), adj


def _forest_ok(adj, S):
    """G[S] acyclic, via |E| == |S| - components."""
    edges = comps = seen = 0
    rest = S
    while rest:
        low = rest & -rest
        stack = [low.bit_length() - 1]
        seen |= low
        comps += 1
        while stack:
            i = stack.pop()
            nb = adj[i] & S
            edges += nb.bit_count()
            nb &= ~seen
            while nb:
                c = nb & -nb
                seen |= c
                stack.append(c.bit_length() - 1)
                nb ^= c
        rest = S & ~seen
    return edges // 2 == S.bit_count() - comps


def _bip_ok(adj, S):
    color = {}
    rest = S
    while rest:
        s = (rest & -rest).bit_length() - 1
        color[s] = 0
        stack = [s]
        while stack:
            i = stack.pop()
            nb = adj[i] & S
            while nb:
                c = nb & -nb
                j = c.bit_length() - 1
                nb ^= c
                if j not in color:
                    color[j] = 1 - color[i]
                    stack.append(j)
                elif color[j] == color[i]:
                    return False
        rest = S & ~sum(1 << k for k in color)
    return True


def _indep_ok(adj, S):
    i, t = 0, S
    while t:
        if t & 1 and (adj[i] & S):
            return False
        t >>= 1
        i += 1
    return True


def _top_down(adj, n, pred):
    for k in range(n, 0, -1):
        for s in itertools.combinations(range(n), k):
            S = 0
            for i in s:
                S |= 1 << i
            if pred(adj, S):
                return k
    return 0


def invariants(g):
    """(f, b, alpha, res) exactly."""
    n, adj = _masks(g)
    return (_top_down(adj, n, _forest_ok),
            _top_down(adj, n, _bip_ok),
            _top_down(adj, n, _indep_ok),
            W.residue(g))


def rhs59(res, b):
    return math.ceil(math.sqrt(res * b))


# ---------------------------------------------------------------------------
# the proved lemmas, as computable predicates
# ---------------------------------------------------------------------------

def max_independent_sets(g):
    n, adj = _masks(g)
    alpha = _top_down(adj, n, _indep_ok)
    return alpha, [S for S in range(1 << n)
                   if S.bit_count() == alpha and _indep_ok(adj, S)]


def lemma_alpha_plus_two(g):
    """L2. Returns alpha+2 when the pair witness exists, else alpha+1.

    PROVED as a lower bound on f: the returned value is the order of an
    exhibited induced forest."""
    n, adj = _masks(g)
    alpha, Is = max_independent_sets(g)
    for I in Is:
        out = [i for i in range(n) if not (I >> i & 1)]
        for x, y in itertools.combinations(out, 2):
            if _forest_ok(adj, I | (1 << x) | (1 << y)):
                return alpha + 2
    return alpha + 1


def toolkit_covers(g, need, with_L2):
    """Does some *structural* (proof-material, non-witness) lower bound on f
    reach `need`?  Exactly wowii_toolkit.FOREST_BOUNDS with the greedy witness
    b_exhibit_forest removed, which is how the published residuals 2 (n<=7) and
    14 (n<=8) were computed -- optionally plus L2.

    Ordered cheapest-first and short-circuiting, because largest_induced_tree is
    exponential and is only reached on the graphs nothing else settles."""
    import wowii_toolkit as T
    n, m = g.number_of_nodes(), g.number_of_edges()
    cheap = [n if nx.is_forest(g) else 0,       # b_acyclic
             nx.diameter(g) + 1,                # b_path
             2 * n - m - 1,                     # b_cyclomatic
             2 * nx.radius(g) - 1]              # b_rad_tree
    if max(cheap) >= need:
        return True
    _, adj = _masks(g)
    alpha = _top_down(adj, n, _indep_ok)
    if alpha + 1 >= need:                       # P2
        return True
    if with_L2 and lemma_alpha_plus_two(g) >= need:
        return True
    b = _top_down(adj, n, _bip_ok)
    if b / 2 + 1 >= need:                       # b_half_bipartite
        return True
    for fn in (T.b_caterpillar_indep, T.b_one_vertex_off, W.largest_induced_tree):
        if fn(g) >= need:
            return True
    return False


def structural_lower_bound(g):
    """A single number for reporting: max over the proved bounds, including L2."""
    n, m = g.number_of_nodes(), g.number_of_edges()
    _, adj = _masks(g)
    alpha = _top_down(adj, n, _indep_ok)
    best = max(alpha + 1, nx.diameter(g) + 1, 2 * n - m - 1,
               2 * nx.radius(g) - 1, lemma_alpha_plus_two(g))
    if nx.is_forest(g):
        best = max(best, n)
    return best


def case_i(f, b, alpha, res):
    """T1: res*b <= (alpha+1)^2."""
    return res * b <= (alpha + 1) ** 2


def case_ii(f, b, alpha, res):
    """T2: b <= f+1."""
    return b <= f + 1


def case_iii(g, f, b, alpha, res):
    """T3: L2 applies and res*b <= (alpha+2)^2."""
    return lemma_alpha_plus_two(g) >= alpha + 2 and res * b <= (alpha + 2) ** 2


def case_i_am(f, b, alpha, res):
    """T1 for the AM-GM form: b + res <= 2*(alpha+1)."""
    return b + res <= 2 * (alpha + 1)


def case_iii_am(g, f, b, alpha, res):
    return lemma_alpha_plus_two(g) >= alpha + 2 and b + res <= 2 * (alpha + 2)


# ---------------------------------------------------------------------------
# the corpus
# ---------------------------------------------------------------------------

def corpus(nmax=8, verbose=True):
    if nmax <= 7:
        return [g for g in nx.graph_atlas_g()
                if 2 <= g.number_of_nodes() <= nmax and nx.is_connected(g)]
    from allgraphs import connected_graphs
    return connected_graphs(nmax, verbose=verbose)


def validate_instrument():
    """The bitmask invariants must agree with wowii.py on every graph the atlas
    holds. Without this the whole file is a report on a different function."""
    bad = 0
    for g in corpus(7):
        f, b, a, r = invariants(g)
        if (f, b, a, r) != (W.largest_induced_forest(g),
                            W.largest_induced_bipartite(g),
                            W.indep_number(g), W.residue(g)):
            bad += 1
            print("   MISMATCH", sorted(g.edges()))
    print(f"   fast invariants vs wowii.py on all 995 graphs n<=7: "
          f"{bad} mismatches" + ("   OK" if not bad else "   *** BROKEN ***"))
    return bad == 0


# ---------------------------------------------------------------------------

def main():
    nmax = 8
    for i, a in enumerate(sys.argv):
        if a == "-n":
            nmax = int(sys.argv[i + 1])

    print(__doc__.split("PROVED")[0].strip().splitlines()[0])
    print()
    print("0. instrument validation")
    if not validate_instrument():
        raise SystemExit("bitmask invariants disagree with wowii.py; stop.")

    print()
    print(f"1. building the corpus (connected graphs on 2..{nmax} vertices)")
    graphs = corpus(nmax)
    print(f"   {len(graphs)} graphs")

    rows = []
    for g in graphs:
        f, b, alpha, res = invariants(g)
        rows.append((g, f, b, alpha, res))

    print()
    print("2. the PROVED lemmas, checked for consistency (a violation would mean"
          " the proof is wrong)")
    checks = {
        "P1  f <= b":                 lambda g, f, b, a, r: f <= b,
        "P2  alpha + 1 <= f":         lambda g, f, b, a, r: a + 1 <= f,
        "P3  b <= 2*alpha":           lambda g, f, b, a, r: b <= 2 * a,
        "P4  f >= ceil(b/2) + 1":     lambda g, f, b, a, r: f >= -(-b // 2) + 1,
        "P5  res <= alpha":           lambda g, f, b, a, r: r <= a,
        "P6  res <= n - ceil(m/Del)": lambda g, f, b, a, r: r <= g.number_of_nodes()
                                      - -(-g.number_of_edges()
                                          // max(d for _, d in g.degree())),
        "L2  alpha+2 witness <= f":   lambda g, f, b, a, r: lemma_alpha_plus_two(g) <= f,
    }
    for name, fn in checks.items():
        bad = sum(1 for row in rows if not fn(*row))
        print(f"   {name:30s} violations {bad:>6}"
              + ("   OK" if not bad else "   *** THE PROOF IS WRONG ***"))

    print()
    print("3. the two statements")
    v59 = [row for row in rows if row[1] < rhs59(row[4], row[2])]
    vam = [row for row in rows if 2 * row[1] < row[2] + row[4]]
    print(f"   C1  conjecture 59      f >= ceil(sqrt(res*b))   violations {len(v59):>6}")
    print(f"   C2  AM-GM strengthening  2*f >= b + res         violations {len(vam):>6}")
    eq = sum(1 for row in rows if 2 * row[1] == row[2] + row[4])
    print(f"       equality 2f = b + res on {eq} of {len(rows)} -- the "
          f"strengthening is sharp")

    print()
    print("4. the case analysis (CHECKED to cover, PROVED case by case)")
    for name, sel in [
            ("(i)  res*b <= (alpha+1)^2", lambda g, f, b, a, r: case_i(f, b, a, r)),
            ("(ii) b <= f+1", lambda g, f, b, a, r: case_ii(f, b, a, r)),
            ("(iii) L2 and res*b <= (alpha+2)^2",
             lambda g, f, b, a, r: case_iii(g, f, b, a, r)),
    ]:
        left = sum(1 for row in rows if not sel(*row))
        print(f"   {name:36s} uncovered {left:>6}")
    left = [row for row in rows
            if not (case_i(*row[1:]) or case_iii(row[0], *row[1:]))]
    print(f"   {'(i) or (iii)':36s} uncovered {len(left):>6}"
          + ("   <-- conjecture 59 settled in range by two cases" if not left else ""))
    left_am = [row for row in rows
               if not (case_i_am(*row[1:]) or case_iii_am(row[0], *row[1:]))]
    print(f"   {'(i) or (iii), AM-GM form':36s} uncovered {len(left_am):>6}"
          + ("   <-- and so is the stronger statement" if not left_am else ""))

    print()
    print("5. the residual under the toolkit's structural bounds "
          "(the published 2 and 14), and what L2 does to it")
    resid = [row for row in rows
             if not toolkit_covers(row[0], rhs59(row[4], row[2]), False)]
    resid_L2 = [row for row in rows
                if not toolkit_covers(row[0], rhs59(row[4], row[2]), True)]
    for cut in sorted({min(7, nmax), nmax}):
        a = sum(1 for row in resid if row[0].number_of_nodes() <= cut)
        c = sum(1 for row in resid_L2 if row[0].number_of_nodes() <= cut)
        print(f"   n <= {cut}: residual without L2 {a:>4},  with L2 {c:>4}")
    print("   the residual graphs, without L2:")
    for g, f, b, alpha, res in resid:
        print(f"      {nx.to_graph6_bytes(g, header=False).strip().decode():>12}"
              f"  n={g.number_of_nodes()} m={g.number_of_edges()}"
              f"  f={f} b={b} alpha={alpha} res={res}"
              f"  rhs={rhs59(res, b)}  girth={0 if nx.is_forest(g) else nx.girth(g)}"
              f"  diam={nx.diameter(g)}"
              f"  degs={sorted((d for _, d in g.degree()), reverse=True)}")
    if resid:
        print(f"   every one of them has res = alpha ("
              f"{sum(1 for r in resid if r[4] == r[3])}/{len(resid)}), "
              f"f = alpha + 2 ({sum(1 for r in resid if r[1] == r[3] + 2)}/{len(resid)}), "
              f"b in {{2*alpha-1, 2*alpha}} "
              f"({sum(1 for r in resid if 2 * r[3] - r[2] in (0, 1))}/{len(resid)}), "
              f"girth 3 ({sum(1 for r in resid if not nx.is_forest(r[0]) and nx.girth(r[0]) == 3)}/{len(resid)})")
        print(f"   and every one holds with EQUALITY f = ceil(sqrt(res*b)) "
              f"({sum(1 for r in resid if r[1] == rhs59(r[4], r[2]))}/{len(resid)})")

    print()
    print("6. C4: where the case analysis provably stops")
    print("   case (iii) needs res*b <= (alpha+2)^2; P3 and P5 allow res*b up to")
    print("   2*alpha^2, and 2*alpha^2 <= (alpha+2)^2 iff alpha <= 4.")
    if resid:
        alphas = sorted({row[3] for row in resid})
        print(f"   the {len(resid)} residual graphs in range have alpha in {alphas}"
              f" -- at the edge")
    wit = nx.from_graph6_bytes(b"JDAp@_PEF??")
    f, b, a, r = invariants(wit)
    print(f"   witness at n=11 (found by hill-climbing): JDAp@_PEF??")
    print(f"      alpha={a} res={r} b={b} f={f}   "
          f"ceil(sqrt(res*b))={rhs59(r, b)}   "
          f"best structural bound={structural_lower_bound(wit)}")
    print(f"      conjecture holds ({f} >= {rhs59(r, b)}); the BOUNDS fail. "
          f"The gap grows.")

    print()
    print("7. C6: tightness of the AM-GM strengthening on K_{a,a}")
    print("      a    f    b  res   2f   b+res")
    for k in range(2, 8):
        g = nx.complete_bipartite_graph(k, k)
        f, b, a_, r = invariants(g)
        print(f"   {k:>4} {f:>4} {b:>4} {r:>4} {2*f:>4} {b+r:>7}"
              f"{'   tight' if 2 * f == b + r else ''}")
    print("   f = a+1 is PROVED (two vertices on each side of K_{a,a} give a")
    print("   4-cycle), b = 2a is immediate, res = 2 is CHECKED below;")
    print("   so the strengthening cannot be improved to 2f >= b + res + 1.")
    bad = [a for a in range(2, 61)
           if W.residue(nx.complete_bipartite_graph(a, a)) != 2]
    print(f"   res(K_a,a) = 2 for a = 2..60: {'confirmed' if not bad else bad}")

    print()
    print("8. FAILED candidates, re-checked so the record cannot rot")
    failed = {
        "F1  f >= alpha + res/2":   lambda g, f, b, a, r: f >= a + r / 2,
        "F2  f >= (alpha + b)/2":   lambda g, f, b, a, r: 2 * f >= a + b,
        "F3  res*b <= (alpha+1)^2": lambda g, f, b, a, r: r * b <= (a + 1) ** 2,
        "F4  b <= 2*res + 2":       lambda g, f, b, a, r: b <= 2 * r + 2,
        "F5  f >= b - res + 1":     lambda g, f, b, a, r: f >= b - r + 1,
    }
    for name, fn in failed.items():
        bad = [row for row in rows if not fn(*row)]
        eg = ""
        if bad:
            sm = min(bad, key=lambda row: row[0].number_of_nodes())
            eg = (f"   smallest n={sm[0].number_of_nodes()} "
                  f"{nx.to_graph6_bytes(sm[0], header=False).strip().decode()}")
        print(f"   {name:28s} violations {len(bad):>6}{eg}")

    print()
    print("Summary: PROVED -- P1..P6, L2, T1..T5, and the two-case analysis is")
    print("valid case by case. CHECKED -- that the two cases cover everything on")
    print("at most 8 vertices, and that 2f >= b + res holds there. NOT PROVED --")
    print("the conjecture. The case analysis is capped at alpha <= 4 and a")
    print("counterexample to the *bounds* exists at n = 11.")


if __name__ == "__main__":
    main()
