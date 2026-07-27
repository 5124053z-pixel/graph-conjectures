r"""Graffiti.pc / Written on the Wall II conjecture 2:  Ls(G) >= 2*(l_avg(G) - 1).

    Ls(G)    = max number of leaves over all spanning trees of G
    l(v)     = alpha(G[N(v)]), the local independence number at v
    l_avg(G) = (1/n) * sum_v l(v)

Everything below separates into three piles, and the piles are kept apart on
purpose.  `main()` re-runs every check.


==============================================================================
PROVED -- holds for every connected graph, by the argument given
==============================================================================

**Master Lemma (leaf count from a forest).**  Let G be connected, n >= 2, and
let F be any acyclic subgraph of G.  Then

    Ls(G)  >=  2 + sum_w max(deg_F(w) - 2, 0).

*Proof.*  Any tree T on n >= 2 vertices has exactly `2 + sum_w (deg_T(w)-2)^+`
leaves: with L leaves and I internal vertices, sum_w deg_T(w) = 2n-2 gives
sum_{internal}(deg_T(w)-2) = 2n-2-L-2I = L-2, and leaves contribute 0 to the
positive part.  F is acyclic and G is connected, so F extends to a spanning
tree T of G, and deg_T(w) >= deg_F(w) for every w.  []

Three corollaries, obtained by choosing F.

**Corollary 1 (classical).**  F = the full star at a maximum-degree vertex
gives `Ls(G) >= Delta(G)`.  This is the bound the project already had.

**THEOREM A (the double-star bound -- the new one).**  For any two distinct
vertices u, v of a connected graph G,

    Ls(G) >= |N(u) u N(v)| - 2,   and  Ls(G) >= |N(u) u N(v)| - 1
                                        when dist(u,v) = 2.

*Proof.*  Three cases, each exhibiting an acyclic F and applying the Master
Lemma.

  uv in E(G).  Put A = N(u)\{v} and B = N(v)\({u} u N(u)).  A and B are
  disjoint, neither contains u or v, and A u B = (N(u) u N(v))\{u,v} -- note
  u in N(v) and v in N(u), so both lie in the union.  Let F consist of the edge
  uv together with ua for a in A and vb for b in B.  F has |A|+|B|+2 vertices
  and |A|+|B|+1 edges and is connected, hence a tree.  deg_F(u) = |A|+1,
  deg_F(v) = |B|+1, so the Master Lemma gives
  Ls >= 2 + (|A|-1)^+ + (|B|-1)^+ >= |A|+|B| = |N(u) u N(v)| - 2.
  (If A or B is empty the inequality only gets slacker.)

  dist(u,v) = 2, with w a common neighbour.  Put A = N(u)\{w},
  B = N(v)\({w} u N(u)).  Now u is not in N(v) and v is not in N(u), so
  A u B = (N(u) u N(v))\{w} and |A|+|B| = |N(u) u N(v)| - 1.  Let F be the
  path u-w-v together with the stars A at u and B at v; again a tree, with
  deg_F(u) = |A|+1, deg_F(v) = |B|+1, deg_F(w) = 2.  Master Lemma gives
  Ls >= |A|+|B| = |N(u) u N(v)| - 1.

  dist(u,v) = d >= 3.  Take a shortest path u = p_0, ..., p_d = v.  Put
  A = N(u)\{p_1}, B = N(v)\{p_{d-1}}.  N(u) n N(v) = empty since d >= 3, so A
  and B are disjoint; no p_j with j >= 2 lies in N(u) and no p_j with
  j <= d-2 lies in N(v).  Let F be the path together with the two stars; a
  tree, deg_F(u) = deg(u), deg_F(v) = deg(v).  Master Lemma gives
  Ls >= deg(u) + deg(v) - 2 = |N(u) u N(v)| - 2.  []

For an edge uv it is convenient to write the bound as
`Ls >= deg(u) + deg(v) - t(uv) - 2`, where t(uv) = |N(u) n N(v)| is the number
of triangles through uv; the two forms are equal because
|N(u) u N(v)| = |N[u] u N[v]| = deg(u) + deg(v) - t(uv) when uv is an edge.

**LEMMA B.**  Let I_v be a maximum independent set of G[N(v)] and let u in I_v.
Then for *any* maximum independent set I_u of G[N(u)],  I_u n I_v = empty, and
therefore  |N(u) u N(v)| >= l(u) + l(v).

*Proof.*  I_v is independent and u in I_v, so no other member of I_v is
adjacent to u, i.e. I_v \ {u} is disjoint from N(u) ⊇ I_u.  And u is not in
I_u since u is not in N(u).  So I_u n I_v = empty.  Both sets lie in
N(u) u N(v).  []

Combining Theorem A and Lemma B:
`Ls(G) >= max_v max_{u in I_v} (l(u) + l(v)) - 2`.

**THEOREM C.  Conjecture 2 is true for every connected triangle-free graph.**

*Proof.*  G triangle-free means N(v) is independent for every v, so
l(v) = deg(v) and l_avg = 2m/n.  Every edge uv lies in no triangle, so
Theorem A reads Ls >= deg(u) + deg(v) - 2.  Averaging over the m edges and
using Cauchy-Schwarz,

    max_{uv in E} (deg(u)+deg(v))  >=  (1/m) sum_{uv in E} (deg u + deg v)
                                    =  (1/m) sum_v deg(v)^2
                                    >= (1/m) (2m)^2 / n  =  4m/n  =  2*l_avg.

Hence Ls >= 2*l_avg - 2 = 2*(l_avg - 1).  []

This is sharp: K_{a,a} has l_avg = a and Ls = 2a - 2, and every cycle C_n with
n >= 4 has l_avg = 2 and Ls = 2.  For K_{a,b} the slack is exactly
(a-b)^2/(a+b).

**THEOREM D (the triangle regime).**  Let G be connected with n >= 2, m >= 1.
If every vertex v that lies on a triangle satisfies deg(v) <= 1 + 2m/n, then
conjecture 2 holds for G.

*Proof.*  Write t(v) for the number of triangles at v, T for the number of
triangles, t(uv) for the number of triangles through the edge uv.  Then
sum_v t(v) = 3T and sum_{uv in E} t(uv) = 3T.

  (i)  By Theorem A,
       Ls + 2 >= max_{uv in E}(deg u + deg v - t(uv))
              >= (1/m) sum_{uv in E}(deg u + deg v - t(uv))
              =  (1/m)(sum_v deg(v)^2 - 3T)
              >= 4m/n - 3T/m                                (Cauchy-Schwarz).

  (ii) G[N(v)] has deg(v) vertices, t(v) edges and maximum degree at most
       deg(v)-1, so any vertex cover of it has size at least t(v)/(deg(v)-1).
       By Gallai, l(v) = alpha(G[N(v)]) = deg(v) - tau(G[N(v)]), so
       l(v) <= deg(v) - t(v)/(deg(v)-1)  (vacuous when t(v) = 0).  Hence
       2*l_avg <= 4m/n - (2/n) sum_v t(v)/(deg(v)-1).

  So it suffices that (2/n) sum_v t(v)/(deg(v)-1) >= 3T/m = (sum_v t(v))/m,
  and that holds term by term as soon as 2/(n(deg(v)-1)) >= 1/m for every v
  with t(v) > 0, i.e. deg(v) <= 1 + 2m/n.  []

Theorem C is the special case T = 0 of Theorem D.

**THE CASE ANALYSIS, and exactly how far it reaches.**

  Case 1.  Delta(G) >= 2*(l_avg - 1).            Then Ls >= Delta settles it.
  Case 2.  every vertex on a triangle has deg(v) <= 1 + 2m/n.   Theorem D.

**PROVED: the two cases are exhaustive whenever 2m/n <= 3**, i.e. for every
connected graph with at most 3n/2 edges.  If both fail then some triangle
vertex has deg(v) > 1 + 2m/n, so Delta > 1 + 2m/n, while Case 1 failing gives
Delta < 2(l_avg - 1) <= 2(2m/n - 1) because l(v) <= deg(v).  Chaining,
1 + 2m/n < 4m/n - 2, i.e. 2m/n > 3.  []

So conjecture 2 is **proved outright** for
  * every triangle-free connected graph, and
  * every connected graph with average degree at most 3, and
  * every connected graph satisfying Case 1 or Case 2 individually.


==============================================================================
VERIFIED COMPUTATIONALLY ONLY -- not proved
==============================================================================

Corpus: all 12,112 connected graphs on 2..8 vertices, generated by
`allgraphs.connected_graphs(8)`, whose counts are checked against OEIS A000088
and A001349 before anything is returned.

  V1.  Conjecture 2 itself: 0 violations over the 12,112.  (Already in the
       repository; re-checked here.)

  V2.  Every bound used above was checked as a lower bound on Ls over all
       12,112 graphs before being relied on: Delta, Theorem A over edges,
       Theorem A over all pairs, and Lemma B's l(u)+l(v)-2.  Zero violations
       each.  Lemma B's own inequality |N(u) u N(v)| >= l(u)+l(v) was checked
       over every (v, I_v, u) triple in the corpus: zero violations.

  V3.  **Case 1 or Case 2 holds for every one of the 12,112 graphs.**  So in
       range the conjecture is settled by proved bounds alone, with no witness
       certificate and no residual.  This is *not* a proof in general: the two
       cases are only provably exhaustive when 2m/n <= 3, and graphs where
       both fail do exist -- see V5.

  V4.  The reduced inequality.  Everything above funnels into one purely local
       statement, with no spanning tree in it:

           (CC)   max_{uv in E(G)}  |N[u] u N[v]|  >=  2 * l_avg(G).

       By Theorem A, CC implies conjecture 2 for every graph.  CC has 0
       violations over the 12,112, and Theorem C's proof is exactly a proof of
       CC for triangle-free graphs.  CC is tight (equality) on C_n for every
       n >= 4, on K_{a,a}, on Q_3, and on the two 3-regular girth-4 graphs on
       8 vertices.  **CC is not proved in general.**

  V5.  The gap regime is non-empty.  Hill-climbing on `min(2(l_avg-1) - Delta,
       max_{v on a triangle} deg(v) - 1 - 2m/n)` finds graphs where both proved
       cases fail, starting at n = 10 (none at n = 9, and none can exist below
       2m/n > 3).  The examples found are stored in `GAP_EXAMPLES`; conjecture 2
       and CC hold on all of them with slack, but they are settled by *nothing
       proved here*.


==============================================================================
THE RESIDUAL
==============================================================================

Before this file, the only structural bound in the project was Ls >= Delta,
whose residual is 32 of the 12,112 (9 of the 995 on <= 7 vertices) -- and the
gap grows with n, which is what made it look like the wrong bound.  Those 32
graphs are characterised as follows, computed in `describe_residual()`:

  * girth is 4 for 28 of them and 3 for the other 4.  **No residual graph has
    girth >= 5**, and none is a tree.
  * 28 of the 32 are triangle-free -- so Theorem C alone removes 28 of them.
  * diameter 2 or 3 for 31 of them (one has diameter 4); radius 2 or 3.
  * near-regular: max degree minus min degree is at most 2 for 30 of them, and
    5 are regular.  Degrees are 3 or 4 with a few 1s and 2s.
  * 17 of the 32 are bipartite.
  * three attain equality Ls = 2(l_avg - 1): K_{3,3}, K_{4,4}, and the cube Q_3.

After Case 1 + Theorem C the residual is **4 graphs**, all on 8 vertices, all
with a triangle: G?rvdo, G?zvfO, G?zveo, GCqjb_.  All four satisfy Case 2, so
after Case 1 + Theorem D the residual is **0** -- in range.


==============================================================================
SEARCHED FOR AND NOT FOUND
==============================================================================

  * `wowii_hunt.py`'s method on conjecture 2 at n = 9..13 (10 restarts x 400
    edge toggles, random G(n,p) starts): best margins -0.11, -0.60, -0.27,
    -0.33, -0.92.  No counterexample.  Note the search never reaches equality
    from random starts, which is why the next item exists.
  * Seeded from the equality family instead (C_n, K_{a,b}, Q_3, Q_4, Petersen,
    Heawood, Franklin), n = 9..11: the margin reaches exactly 0 -- the
    unperturbed cycles -- and never exceeds it.  **This corrects a claim in the
    README, that conjecture 2 "never reaches equality" at these orders.  It
    does: on every C_n with n >= 4, on K_{a,a}, and on the cube Q_3.**  The
    earlier reading came from measuring only random-start hill climbs, which
    do not go near the tight graphs.  Rerun with `--hunt 9-13`.
  * CC hunted directly (40 restarts x 1500 toggles) at n = 9..18: best margins
    -0.11, 0.00, 0.00, 0.00, -0.08, 0.00, -0.07, 0.00, -0.24, -0.22.  Equality
    is reached at n = 10, 11, 12, 14 and 16 and never exceeded.  No
    counterexample to CC.  CC is the cheapest thing here to refute -- it needs
    no spanning tree, only neighbourhood independence numbers -- so this is the
    strongest of the three searches, and it still found nothing.

How much that is worth: conjecture 200 of this same list looked settled over
all 12,112 and is false at 11 vertices, and this project's own hunt failed to
rediscover that counterexample.  A null hunt is weak evidence.  What is *not*
weak is Theorems A, C and D, which hold for all graphs.


==============================================================================
WHAT DID NOT WORK -- recorded rather than hidden
==============================================================================

  * Proving CC in general.  Averaging over edges with uniform weights, then
    Cauchy-Schwarz, is what proves the triangle-free case; with triangles the
    same route needs `sum_v l(v)*c(v) >= (sum l)(sum c)/n` for c the in-degree
    of the arc set {v -> u : u in I_v}, which is Chebyshev's sum inequality and
    needs l and c similarly ordered.  Nothing forces that.  Theorem D is the
    salvage: it buys the triangle correction with the crude vertex-cover
    estimate tau >= e/Delta, and pays for it with the hypothesis
    deg(v) <= 1 + 2m/n.
  * Lemma B's bound `l(u)+l(v)-2` on its own is *weaker* than Delta on dense
    graphs: it leaves a residual of 28 graphs (all girth 3, all with Delta
    large), so it is only useful in combination.  It is not a replacement for
    Ls >= Delta.
  * Reducing to a spanning triangle-free subgraph.  Ls is monotone under edge
    addition (a connected dominating set survives adding edges), so it would
    suffice to find a connected spanning triangle-free H with l_H >= l_G
    pointwise, and then apply Theorem C to H.  No such construction was found,
    and it was not tested exhaustively -- searching spanning subgraphs is 2^m.
  * `Ls >= Delta + 1` on the residual is a real pattern (all 9 of the <= 7
    residual graphs have Ls = Delta + 1) and is false in general -- stars and
    cycles have Ls = Delta.  It was not pursued.
"""
from __future__ import annotations

import argparse
import itertools
import math
import time

import networkx as nx

import wowii as W


# ---------------------------------------------------------------------------
# invariants
# ---------------------------------------------------------------------------

def Ls(g):
    """Maximum leaf number, via n - (minimum connected dominating set)."""
    return W.max_leaves_via_cds(g)


def lvals(g):
    return {v: W.local_independence(g, v) for v in g}


def l_avg(g):
    return W.avg_l(g)


def rhs(g):
    """The conjectured lower bound 2*(l_avg - 1)."""
    return 2 * (l_avg(g) - 1)


def girth(g):
    return 0 if nx.is_forest(g) else nx.girth(g)


def max_independent_sets(h):
    """Every maximum independent set of h, as frozensets."""
    if h.number_of_nodes() == 0:
        return []
    cliques = list(nx.find_cliques(nx.complement(h)))
    k = max(len(c) for c in cliques)
    return [frozenset(c) for c in cliques if len(c) == k]


# ---------------------------------------------------------------------------
# the proved lower bounds on Ls
# ---------------------------------------------------------------------------

def b_master(g, F):
    """Master Lemma: 2 + sum_w (deg_F(w) - 2)^+ for an acyclic subgraph F.

    F is given as an iterable of edges of g. Raises if F has a cycle, because
    a silently-wrong forest would turn this into an unsound bound."""
    h = nx.Graph()
    h.add_nodes_from(g)
    h.add_edges_from(F)
    if not nx.is_forest(h):
        raise ValueError("b_master needs an acyclic F")
    return 2 + sum(max(d - 2, 0) for _, d in h.degree())


def b_maxdeg(g):
    """Corollary 1 (classical): Ls >= Delta. The Master Lemma with F the full
    star at a maximum-degree vertex."""
    return max(d for _, d in g.degree())


def b_edge_ball(g):
    """Theorem A over edges: Ls >= |N[u] u N[v]| - 2."""
    return max(len(set(g[u]) | set(g[v]) | {u, v}) for u, v in g.edges())- 2


def b_pair_ball(g):
    """Theorem A over all pairs, with the distance-2 case worth one more."""
    dist = dict(nx.all_pairs_shortest_path_length(g))
    best = 0
    for u, v in itertools.combinations(list(g), 2):
        union = len(set(g[u]) | set(g[v]))
        best = max(best, union - (1 if dist[u][v] == 2 else 2))
    return best


def b_lpair(g, L=None):
    """Theorem A + Lemma B: Ls >= max_v max_{u in I_v} (l(u) + l(v)) - 2."""
    L = L or lvals(g)
    best = 0
    for v in g:
        nb = list(g[v])
        if not nb:
            continue
        for I in max_independent_sets(g.subgraph(nb)):
            for u in I:
                best = max(best, L[u] + L[v])
    return best - 2


PROVED_BOUNDS = [b_maxdeg, b_edge_ball, b_pair_ball, b_lpair]


# ---------------------------------------------------------------------------
# the two proved regimes
# ---------------------------------------------------------------------------

def case1(g):
    """Delta >= 2*(l_avg - 1).  Settled by Ls >= Delta."""
    return b_maxdeg(g) >= rhs(g) - 1e-9


def case2(g):
    """Theorem D: every vertex on a triangle has deg(v) <= 1 + 2m/n."""
    n, m = g.number_of_nodes(), g.number_of_edges()
    if m == 0:
        return True
    tri = nx.triangles(g)
    return all(g.degree(v) <= 1 + 2 * m / n + 1e-12 for v in g if tri[v] > 0)


def cc(g):
    """The reduced local inequality: max_e |N[u] u N[v]| >= 2*l_avg.

    VERIFIED IN RANGE ONLY.  Proved for triangle-free graphs (Theorem C) and
    under Theorem D's hypothesis; open in general."""
    return b_edge_ball(g) + 2 >= 2 * l_avg(g) - 1e-9


# ---------------------------------------------------------------------------
# graphs found in the regime the case analysis does not reach (V5)
# ---------------------------------------------------------------------------

#: Both Case 1 and Case 2 fail on these.  Found by hill-climbing; none exists
#: on <= 9 vertices.  Conjecture 2 holds on all of them, by direct computation
#: and by nothing proved in this file.
GAP_EXAMPLES = [
    "IWGe}zLv?",
    "IreSdHfZo",
    "I]VE@{{N_",
    "JlXBD[]h_{?",
    "JrG[rKXeJL?",
    "K[aXIvEU`oCm",
    "L{D|[XZhJunYYs",
]


# ---------------------------------------------------------------------------
# the corpus
# ---------------------------------------------------------------------------

def corpus(nmax=8, verbose=False):
    if nmax <= 7:
        return [g for g in nx.graph_atlas_g()
                if 2 <= g.number_of_nodes() <= nmax and nx.is_connected(g)]
    from allgraphs import connected_graphs
    return connected_graphs(nmax, verbose=verbose)


# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------

def check_instrument(graphs):
    """Ls via connected dominating sets must agree with brute-force spanning
    tree enumeration.  Only run on <= 6 vertices: K_8 alone has 8^6 spanning
    trees."""
    bad = [g for g in graphs if g.number_of_nodes() <= 6
           and Ls(g) != W.max_leaves_spanning_tree(g)]
    small = sum(1 for g in graphs if g.number_of_nodes() <= 6)
    return small, bad


def check_master_lemma(graphs):
    """The Master Lemma, checked on an explicitly built forest per graph: the
    double star at a maximum-degree vertex and one of its neighbours."""
    bad = []
    for g in graphs:
        v = max(g, key=g.degree)
        u = next(iter(g[v]))
        A = set(g[v]) - {u}
        B = set(g[u]) - {v} - A
        F = [(u, v)] + [(v, a) for a in A] + [(u, b) for b in B]
        if b_master(g, F) > Ls(g):
            bad.append(g)
    return bad


def check_bounds(graphs):
    """Every bound must be <= Ls on every graph.  A bound that fails once is
    unsound and is discarded, not patched."""
    out = {}
    for b in PROVED_BOUNDS:
        bad = [g for g in graphs if b(g) > Ls(g)]
        out[b.__name__] = bad
    return out


def check_lemma_b(graphs):
    """|N(u) u N(v)| >= l(u) + l(v) for every v, every maximum independent set
    I_v of G[N(v)], and every u in I_v."""
    bad = []
    for g in graphs:
        L = lvals(g)
        for v in g:
            nb = list(g[v])
            if not nb:
                continue
            for I in max_independent_sets(g.subgraph(nb)):
                for u in I:
                    if len(set(g[u]) | set(g[v])) < L[u] + L[v]:
                        bad.append((g, u, v))
    return bad


def check_theorem_c(graphs):
    """On triangle-free graphs: l(v) == deg(v) for every v, and the
    Cauchy-Schwarz step max_e (deg u + deg v) >= 4m/n."""
    bad_l, bad_cs = [], []
    tf = [g for g in graphs if sum(nx.triangles(g).values()) == 0]
    for g in tf:
        L = lvals(g)
        if any(L[v] != g.degree(v) for v in g):
            bad_l.append(g)
        n, m = g.number_of_nodes(), g.number_of_edges()
        if max(g.degree(u) + g.degree(v) for u, v in g.edges()) < 4 * m / n - 1e-9:
            bad_cs.append(g)
    return len(tf), bad_l, bad_cs


def check_case_analysis(graphs):
    left = [g for g in graphs if not (case1(g) or case2(g))]
    return left


def check_conjecture(graphs):
    return [g for g in graphs if Ls(g) < rhs(g) - 1e-9]


def check_cc(graphs):
    return [g for g in graphs if not cc(g)]


def describe_residual(graphs):
    """The 32 graphs where Ls >= Delta alone does not reach 2*(l_avg - 1)."""
    res = [g for g in graphs if not case1(g)]
    rows = []
    for g in res:
        degs = sorted(d for _, d in g.degree())
        rows.append(dict(
            g6=nx.to_graph6_bytes(g, header=False).strip().decode(),
            n=g.number_of_nodes(), m=g.number_of_edges(), degs=degs,
            girth=girth(g), diam=nx.diameter(g), rad=nx.radius(g),
            bip=nx.is_bipartite(g),
            tf=sum(nx.triangles(g).values()) == 0,
            regular=degs[0] == degs[-1], spread=degs[-1] - degs[0],
            Ls=Ls(g), rhs=rhs(g), lavg=l_avg(g),
        ))
    return res, rows


def check_tight_family():
    """Equality cases.  Ls = 2*(l_avg - 1) exactly, so the conjecture cannot be
    improved by any additive constant."""
    out = []
    for n in range(4, 13):
        g = nx.cycle_graph(n)
        out.append((f"C_{n}", Ls(g), rhs(g)))
    for a in range(2, 6):
        g = nx.complete_bipartite_graph(a, a)
        out.append((f"K_{a},{a}", Ls(g), rhs(g)))
    q3 = nx.convert_node_labels_to_integers(nx.hypercube_graph(3))
    out.append(("Q_3", Ls(q3), rhs(q3)))
    return out


def check_kab_slack():
    """K_{a,b} with 2 <= a <= b has slack exactly (a-b)^2/(a+b).

    l = b on the a-side and a on the b-side, so l_avg = 2ab/(a+b); two adjacent
    vertices dominate everything so gamma_c = 2 and Ls = a+b-2; the slack is
    ((a+b)^2 - 4ab)/(a+b).  a = 1 is excluded on purpose: the star has a
    dominating vertex, gamma_c = 1 and Ls = b, and the formula does not apply.
    That exception was found by this check failing, not by assuming."""
    bad = []
    for a in range(2, 7):
        for b in range(a, 7):
            g = nx.complete_bipartite_graph(a, b)
            got = Ls(g) - rhs(g)
            want = (a - b) ** 2 / (a + b)
            if abs(got - want) > 1e-9:
                bad.append((a, b, got, want))
    return bad


def check_gap_examples():
    """The graphs where BOTH proved cases fail.  The conjecture and CC are
    checked directly on them; nothing in this file proves them."""
    rows = []
    for s in GAP_EXAMPLES:
        g = nx.from_graph6_bytes(s.encode())
        rows.append(dict(g6=s, n=g.number_of_nodes(), m=g.number_of_edges(),
                         case1=case1(g), case2=case2(g),
                         davg=2 * g.number_of_edges() / g.number_of_nodes(),
                         holds=Ls(g) >= rhs(g) - 1e-9, cc=cc(g)))
    return rows


# ---------------------------------------------------------------------------
# the searches whose null results the docstring reports.  Kept here so the
# claims can be re-run rather than taken on trust.
# ---------------------------------------------------------------------------

def _conj_margin(g):
    """> 0  <=>  counterexample to conjecture 2."""
    if not nx.is_connected(g) or g.number_of_nodes() < 2:
        return None
    return rhs(g) - Ls(g)


def _cc_margin(g):
    """> 0  <=>  counterexample to CC."""
    if not nx.is_connected(g) or g.number_of_edges() == 0:
        return None
    return 2 * l_avg(g) - (b_edge_ball(g) + 2)


def _climb(g, margin_fn, iters, rng):
    cur = margin_fn(g)
    if cur is None:
        return g, -99.0
    best = (g, cur)
    for _ in range(iters):
        u, v = rng.sample(list(g), 2)
        h = g.copy()
        h.remove_edge(u, v) if h.has_edge(u, v) else h.add_edge(u, v)
        m = margin_fn(h)
        if m is None:
            continue
        if m >= cur:
            g, cur = h, m
            if cur > best[1]:
                best = (g, cur)
        if cur > 1e-9:
            return g, cur
    return best


def tight_seeds(n):
    """The equality family, and the near-misses around it.  Random G(n,p)
    starts never climb to equality on conjecture 2 -- every reported margin in
    wowii_hunt is strictly negative -- so the only sensible starting points are
    the graphs that already attain it."""
    out = []
    if n >= 4:
        out.append((f"C_{n}", nx.cycle_graph(n)))
    for a in range(2, n // 2 + 1):
        out.append((f"K_{a},{n-a}", nx.complete_bipartite_graph(a, n - a)))
    if n in (8, 16):
        out.append((f"Q_{3 if n == 8 else 4}",
                    nx.hypercube_graph(3 if n == 8 else 4)))
    if n == 10:
        out.append(("Petersen", nx.petersen_graph()))
    if n == 14:
        out.append(("Heawood", nx.heawood_graph()))
    if n == 12:
        out.append(("Franklin", nx.LCF_graph(12, [5, -5], 6)))
    return [(nm, nx.convert_node_labels_to_integers(g)) for nm, g in out]


def hunt(lo=9, hi=13, which="all"):
    """Report only.  A null hunt is weak evidence: conjecture 200 of this list
    looked settled over all 12,112 graphs and is false at 11 vertices."""
    import random
    for n in range(lo, hi + 1):
        line = [f"  n={n:>3}"]
        if which in ("all", "random"):
            best = -99.0
            for seed in range(10):
                rng = random.Random(seed * 7919 + n)
                p = rng.choice([0.2, 0.3, 0.4, 0.5, 0.6])
                g = nx.gnp_random_graph(n, p, seed=rng.randrange(10 ** 6))
                if not nx.is_connected(g):
                    continue
                _, m = _climb(g, _conj_margin, 400, rng)
                best = max(best, m)
                if m > 1e-9:
                    break
            line.append(f"conj/random {best:+.4f}")
        if which in ("all", "tight"):
            best, who = -99.0, ""
            for nm, s in tight_seeds(n):
                m0 = _conj_margin(s)
                if m0 is not None and m0 > best:
                    best, who = m0, nm
                for seed in range(4):
                    rng = random.Random(seed * 104729 + n)
                    _, m = _climb(s.copy(), _conj_margin, 200, rng)
                    if m > best:
                        best, who = m, nm
            line.append(f"conj/tight {best:+.4f} ({who})")
        if which in ("all", "cc"):
            best = -99.0
            for seed in range(15):
                rng = random.Random(seed * 7919 + n)
                p = rng.choice([0.15, 0.3, 0.4, 0.5, 0.6, 0.7])
                g = nx.gnp_random_graph(n, p, seed=rng.randrange(10 ** 6))
                if not nx.is_connected(g) or g.number_of_edges() == 0:
                    continue
                _, m = _climb(g, _cc_margin, 800, rng)
                best = max(best, m)
                if m > 1e-9:
                    break
            line.append(f"CC {best:+.4f}")
        flag = "   *** COUNTEREXAMPLE ***" if any(
            float(x.split()[1]) > 1e-9 for x in line[1:]) else ""
        print("   ".join(line) + flag, flush=True)


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=8,
                    help="scan all connected graphs on 2..n vertices (default 8)")
    ap.add_argument("--hunt", nargs="?", const="9-13", default=None,
                    metavar="LO-HI",
                    help="also hunt above the exhaustive range, e.g. --hunt 9-13")
    args = ap.parse_args()

    t0 = time.time()
    print(f"generating all connected graphs on 2..{args.n} vertices")
    graphs = corpus(args.n, verbose=args.n > 7)
    print(f"  {len(graphs)} graphs   [{time.time() - t0:.0f}s]\n")

    ok = True

    def report(label, bad, extra=""):
        nonlocal ok
        if bad:
            ok = False
        print(f"  {'OK ' if not bad else '*** FAIL ***'} {label:<58}"
              f"{'0 violations' if not bad else str(len(bad)) + ' VIOLATIONS'}{extra}")

    print("instrument")
    small, bad = check_instrument(graphs)
    report(f"Ls via CDS == brute force, on the {small} graphs with n <= 6", bad)

    print("\nproved bounds -- every one must be <= Ls on every graph")
    report("Master Lemma on an explicit double-star forest",
           check_master_lemma(graphs))
    for name, bad in check_bounds(graphs).items():
        report(f"{name}", bad)
    report("Lemma B: |N(u) u N(v)| >= l(u)+l(v) for u in I_v",
           check_lemma_b(graphs))

    print("\nTheorem C (triangle-free)")
    ntf, bad_l, bad_cs = check_theorem_c(graphs)
    report(f"l(v) == deg(v) on the {ntf} triangle-free graphs", bad_l)
    report("max_e (deg u + deg v) >= 4m/n  (Cauchy-Schwarz step)", bad_cs)

    print("\nthe conjecture, and the case analysis")
    report("conjecture 2 itself", check_conjecture(graphs))
    left = check_case_analysis(graphs)
    report("every graph satisfies Case 1 or Case 2 (proved regimes)", left)
    report("reduced inequality CC: max_e |N[u] u N[v]| >= 2*l_avg",
           check_cc(graphs))

    print("\nresidual of Ls >= Delta alone (what this file had to close)")
    res, rows = describe_residual(graphs)
    print(f"    {len(res)} graphs, by order: "
          f"{ {k: sum(1 for r in rows if r['n'] == k) for k in sorted({r['n'] for r in rows})} }")
    print(f"    girth: { {k: sum(1 for r in rows if r['girth'] == k) for k in sorted({r['girth'] for r in rows})} }"
          f"   -- none has girth >= 5")
    print(f"    triangle-free: {sum(1 for r in rows if r['tf'])} of {len(rows)}"
          f"  (Theorem C removes exactly these)")
    print(f"    bipartite: {sum(1 for r in rows if r['bip'])}   "
          f"regular: {sum(1 for r in rows if r['regular'])}   "
          f"degree spread <= 2: {sum(1 for r in rows if r['spread'] <= 2)}")
    print(f"    diameter: { {k: sum(1 for r in rows if r['diam'] == k) for k in sorted({r['diam'] for r in rows})} }"
          f"   radius: { {k: sum(1 for r in rows if r['rad'] == k) for k in sorted({r['rad'] for r in rows})} }")
    eq = [r for r in rows if abs(r["Ls"] - r["rhs"]) < 1e-9]
    print(f"    at equality Ls = 2(l_avg - 1): {[r['g6'] for r in eq]}")
    after_c = [r for r in rows if not r["tf"]]
    print(f"    residual after Case 1 + Theorem C: {len(after_c)} "
          f"{[r['g6'] for r in after_c]}")
    still = [g for g in res if not case2(g)]
    report("residual after Case 1 + Theorem D", still)

    print("\ntightness -- the conjecture is sharp, so no additive slack exists")
    tight = check_tight_family()
    bad = [(nm, a, b) for nm, a, b in tight if abs(a - b) > 1e-9]
    report(f"Ls == 2(l_avg - 1) exactly on {len(tight)} graphs "
           f"(C_4..C_12, K_a,a, Q_3)", bad)
    report("K_{a,b} slack is exactly (a-b)^2/(a+b)", check_kab_slack())

    print("\nthe regime the case analysis does NOT reach (V5)")
    for r in check_gap_examples():
        flag = "" if (r["holds"] and r["cc"]) else "   *** LOOK HERE ***"
        if not (r["holds"] and r["cc"]):
            ok = False
        print(f"    {r['g6']:>16} n={r['n']:>2} 2m/n={r['davg']:.2f} "
              f"case1={r['case1']} case2={r['case2']} "
              f"conjecture holds={r['holds']} CC holds={r['cc']}{flag}")
    print("    these are settled by direct computation only; the proofs above")
    print("    say nothing about them.  None exists below 10 vertices.")

    if args.hunt:
        lo, hi = (int(x) for x in args.hunt.split("-"))
        print(f"\ncounterexample hunt, {lo}..{hi} vertices "
              f"(positive margin = counterexample)")
        hunt(lo, hi)

    print(f"\n{'all checks passed' if ok else '*** SOMETHING FAILED ***'}"
          f"   [{time.time() - t0:.0f}s]")
    print("\nPROVED for all graphs : Master Lemma, Theorem A, Lemma B,")
    print("                        Theorem C (triangle-free), Theorem D,")
    print("                        and Case 1 u Case 2 exhaustive when 2m/n <= 3.")
    print("VERIFIED in range only: Case 1 u Case 2 covering everything, and CC.")
    print("NOT PROVED             : the conjecture in general.")


if __name__ == "__main__":
    main()
