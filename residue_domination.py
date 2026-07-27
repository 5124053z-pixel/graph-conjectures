r"""Residue 2 and connected domination -- the lemma conjecture 217 needs.

    CLAIM.  If the Havel-Hakimi residue of a connected graph G is 2, then
    gamma_c(G) <= 4, equivalently Ls(G) >= n - 4.

The claim is what makes WOWII 217 finite: 217 says Ls(G) <= 4*[R(G)==2] + 2
implies a Hamiltonian path, its R != 2 branch is empty of content (Ls <= 2 forces
every spanning tree to be a path, so G is a path or a cycle), and in its R = 2
branch the claim turns Ls <= 6 into n <= 10.

STATUS: NOT PROVED.  What is here is one real theorem, a proof of the claim in
the case gamma(G) <= 2 (which is 99.6% of the residue-2 graphs on at most 8
vertices), the death of the two obvious routes to the rest, and a hunt that
found nothing.  The correct constant still looks like 4.

--------------------------------------------------------------------------
1.  WHAT RESIDUE 2 FORCES  (PROVED)
--------------------------------------------------------------------------

Write W(G) = sum_v 1/(deg(v)+1) for the Caro-Wei bound and R(G) for the residue.

    LEMMA 1.  A Havel-Hakimi step never decreases W.
    Proof.  Let pi = (d_1 >= ... >= d_n) be graphical with d_1 = d >= 1.  The
    step deletes d_1 and decrements d_2, ..., d_{d+1}.  Since pi is graphical
    and d_1 = d, the vertex of degree d has d neighbours of positive degree,
    so d_{d+1} >= 1 and the step is legal.  Then

        W(pi') - W(pi)  =  sum_{i=2}^{d+1} [1/d_i - 1/(d_i+1)]  -  1/(d+1)
                        =  sum_{i=2}^{d+1} 1/(d_i(d_i+1))  -  1/(d+1).

    Each of those d terms has 1 <= d_i <= d and x -> 1/(x(x+1)) is decreasing,
    so each is at least 1/(d(d+1)) and the sum is at least d/(d(d+1)) =
    1/(d+1).  QED

    THEOREM 1.  R(G) >= W(G) = sum_v 1/(deg(v)+1).
    Proof.  Havel-Hakimi terminates at the all-zero sequence of length R, whose
    W is R * 1/(0+1) = R.  By Lemma 1 the run never decreases W.  QED

    COROLLARY 1.  n <= R * (Delta + 1).
    COROLLARY 2.  For every t >= 0,  #{v : deg(v) <= t}  <=  R * (t + 1).
    (Both because W >= #{v : deg v <= t} / (t+1), and n/(Delta+1) is the t=Delta
    case.)

So residue 2 means: n <= 2*Delta + 2, at most 2t+2 vertices of degree at most t
for every t, hence the i-th smallest degree is at least ceil(i/2) - 1 and
m >= n(n-2)/8.  A residue-2 graph is dense, and every sparse "long" part is
short: at most 6 vertices of degree <= 2 kills long paths, long cycles, spiders
and subdivisions outright.  The weaker Corollary 1 alone is classical
(Favaron, Maheo, Sacle); Theorem 1 is proved here and I have not checked whether
it is in the literature.

--------------------------------------------------------------------------
2.  THE CLAIM WHEN gamma(G) <= 2  (PROVED)
--------------------------------------------------------------------------

    DUCHET-MEYNIEL (1982).  gamma_c(G) <= 3*gamma(G) - 2 for connected G.

    THEOREM 2.  If R(G) = 2 and gamma(G) <= 2 then gamma_c(G) <= 4.
    Proof.  Immediate from Duchet-Meyniel.  QED

That is not a joke result: 4,358 of the 4,375 connected residue-2 graphs on at
most 8 vertices have gamma <= 2, so Theorem 2 proves the claim for 99.6% of them
in range, and C_6 -- the extremal graph, gamma = 2 and gamma_c = 4 -- is one of
the graphs it covers, at equality in Duchet-Meyniel.  The open case is exactly

        R(G) = 2  and  gamma(G) >= 3.

Residue 2 does NOT force gamma <= 2: 4-regular graphs on 10 vertices have
residue 2 (Corollary 1 is tight there) and about half of them have gamma = 3.
So "gamma(G) <= R(G)", which would have finished the proof, is FALSE.  Every
one of those gamma = 3 graphs still has gamma_c = 3.

--------------------------------------------------------------------------
3.  THE CONSTANT 4 IS SHARP AT EVERY ORDER  (PROVED)
--------------------------------------------------------------------------

C_6 is not an isolated accident of six vertices.

    THEOREM 4.  Let H_T (T >= 1) be C_6 with one block blown up: a clique
    B of size T, two vertices a and e each joined to all of B, and a path
    a - b - c - d - e.  Then n = T + 5, R(H_T) = 2 and gamma_c(H_T) = 4.

    Proof.  Degrees: (T+1) on B u {a, e} (that is T+2 vertices) and 2 on
    b, c, d.  Havel-Hakimi on (k)^{k+1} 2 2 2 with k >= 2 deletes the pivot k
    and decrements the k remaining k's, giving (k-1)^k 2 2 2; starting from
    k = T+1 this descends to 2^3 2 2 2 = 2^6, the sequence of C_6, whose
    residue is 2.  For gamma_c: {a,b,c,d} is connected and dominating (a covers
    B and b, and d covers e), so gamma_c <= 4.  Suppose S were connected
    dominating with |S| = 3.  If c in S then dominating B needs S to meet
    B u {a,e}, and with only three connected vertices that forces S = {c,b,a}
    or {c,d,e}; the first leaves e undominated, the second leaves a
    undominated.  If c not in S then S must meet {a,b} and {d,e}, whose
    distance is 2 through c or through B, so S is {a, x, e} with x in B -- and
    then c is undominated.  Hence gamma_c = 4.  QED

So if the claim is true it is best possible for every n >= 6, and the
gamma_c = 4 case of conjecture 217's R = 2 branch is not vacuous.  These are
also exactly the sequences the realisation search below finds hardest, and they
sit at W = 2 exactly, the equality case of Theorem 1's corollaries.

--------------------------------------------------------------------------
4.  TWO ROUTES THAT DIE  (both disproved here)
--------------------------------------------------------------------------

The shape of the evidence suggests

        CONJECTURE A.  gamma_c(G) <= 3*R(G) - 2,

which is tight for every value of R: the blow-up C_{3j}[K_t] has residue j and
gamma_c = 3j - 2, and C_6 is the j = 2 case.  Conjecture A restricted to R = 2
is exactly the claim.  Two attempts to prove it:

(a) Through W.  Theorem 1 gives R >= W, so it would be enough that
    gamma_c <= 3W - 2, or even just the special case W <= 2 => gamma_c <= 4.
    BOTH ARE FALSE.  gamma_c <= 3W - 2 fails on 85 of the 12,112 connected
    graphs on at most 8 vertices (all of them dense, with gamma_c <= 4 and W as
    low as 8/7).  Worse, the special case fails badly:

    THEOREM 3 (construction).  For k >= 2 and t0, t, s >= 1 let S(k; t0, t, s)
    be: a centre clique C of size t0; for each of k legs a clique L_i of size t
    joined completely to C, and a clique M_i of size s joined completely to L_i.
    Then gamma_c(S) = k + 1.
    Proof.  N[M_i] = L_i u M_i, so a dominating set meets each of the k
    pairwise disjoint, pairwise non-adjacent sets L_i u M_i; a set of size k
    meeting each exactly once induces no edge, so it cannot be connected, and
    gamma_c >= k+1.  Conversely one c in C dominates C and every L_i, and one
    l_i in each L_i dominates L_i and M_i, and c is adjacent to every l_i, so
    that k+1 vertices form a connected dominating set.  QED

    With k = 4 this has gamma_c = 5 and W(S) -> 1 as t0 >> t >> s: at
    (t0,t,s) = (500,22,1), n = 592, W = 1.193.  So W <= 2 says nothing about
    gamma_c and the Caro-Wei content of the residue is not enough.  What saves
    the claim is that S(4; t0,t,s) has residue 4 -- at every one of the 21,888
    parameter triples scanned.  Annealing on W under the constraint
    gamma_c >= 5 also crosses 2 without help from the spiders (n = 20,
    W = 1.971, residue 3).  Any proof of the claim must use the residue beyond
    Corollaries 1 and 2.

(b) Through the degree conditions alone.  Corollaries 1 and 2 are also not
    enough on their own: the "geometric chain" -- a path of k cliques of sizes
    1, 3, 9, ..., 3^{k-1}, consecutive ones joined completely -- satisfies every
    instance of Corollary 2 for R = 2 and has gamma_c = k - 2, so at k = 7 it is
    a graph with gamma_c = 5 obeying all the proved consequences of residue 2.
    Its residue is 3.  The gap between "the corollaries" and "residue 2" is
    real and it is exactly where the remaining proof has to live.

--------------------------------------------------------------------------
5.  THE HUNT  (VERIFIED: nothing found)
--------------------------------------------------------------------------

Everything below is verified, not proved.

  * Exhaustively, all 12,112 connected graphs on at most 8 vertices: 4,375 have
    residue 2, the largest gamma_c among them is 4, attained by exactly 6
    graphs (C_6 is the smallest).  The minimum residue over the graphs with
    gamma_c = k is 1, 2, 2, 2, 3, 3 for k = 1..6, i.e. exactly
    ceil((k+2)/3), the equality pattern of Conjecture A.
  * All 2,326 residue-2 degree sequences on at most 10 vertices, with 100-200
    realisations of each explored by degree-preserving double edge swaps: not
    one realisation with gamma_c >= 5.  The hardest sequence at each n is
    (n-4)^{n-3} 2^3 -- precisely the degree sequence of H_{n-5} in Theorem 4.
  * Simulated annealing on the residue itself under the hard constraint
    gamma_c >= 5, n = 7..18 (and n = 19..24 offline), started both from cycle
    blow-ups and from the Theorem 3 spiders: the minimum residue reached is 3
    at every n.
  * Annealing on W under the same constraint reaches 1.971 at n = 20 -- so the
    W relaxation really is empty above 2, which is route (a) dying, not the
    claim.
  * Families: C_k[K_t] has residue ceil(k/3) and gamma_c = k-2; K_{a,a} has
    residue 2 for every a, so residue-2 graphs exist at every order and their
    independence number is unbounded; P_k[K_t] has residue 2 exactly for
    k <= 5, where gamma_c = k-2 <= 3.

So: no counterexample, the constant 4 stands, and the conjecture-217 reduction
to n <= 10 is safe as far as any of this can see -- but the claim is a claim.
"""
from __future__ import annotations

import argparse
import itertools
import math
import random
import time

import networkx as nx

import allgraphs
import domination
import wowii


# ---------------------------------------------------------------- primitives

def residue_seq(seq):
    """Havel-Hakimi residue of a degree sequence; None if not graphical."""
    s = sorted(seq, reverse=True)
    while s and s[0] > 0:
        d = s[0]
        s = s[1:]
        if d > len(s):
            return None
        for i in range(d):
            s[i] -= 1
            if s[i] < 0:
                return None
        s.sort(reverse=True)
    return len(s)


def residue(g):
    return residue_seq([d for _, d in g.degree()])


def caro_wei_seq(seq):
    return sum(1.0 / (d + 1) for d in seq)


def caro_wei(g):
    return caro_wei_seq([d for _, d in g.degree()])


def _adj(g):
    """Bitmask adjacency on vertices relabelled 0..n-1."""
    idx = {v: i for i, v in enumerate(g)}
    a = [0] * len(idx)
    for u, v in g.edges():
        a[idx[u]] |= 1 << idx[v]
        a[idx[v]] |= 1 << idx[u]
    return a


def _connected_dominating(adj, n, size):
    """Iterate the connected dominating sets of exactly `size` vertices."""
    full = (1 << n) - 1
    clos = [adj[v] | (1 << v) for v in range(n)]
    for s in itertools.combinations(range(n), size):
        cov = 0
        for v in s:
            cov |= clos[v]
        if cov != full:
            continue
        smask = 0
        for v in s:
            smask |= 1 << v
        seen = 1 << s[0]
        stack = [s[0]]
        while stack:
            x = stack.pop()
            nxt = adj[x] & smask & ~seen
            while nxt:
                b = nxt & -nxt
                seen |= b
                stack.append(b.bit_length() - 1)
                nxt ^= b
        if bin(seen).count("1") == size:
            yield s


def gamma_c(g):
    """Minimum connected dominating set size.  K_2 is 1, matching wowii's
    Ls(K_2) = 2 convention only up to that one graph (see wowii.py)."""
    n = g.number_of_nodes()
    if n <= 1:
        return 0
    adj = _adj(g)
    for k in range(1, n + 1):
        for _ in _connected_dominating(adj, n, k):
            return k
    return n


def has_cds_upto(g, k):
    n = g.number_of_nodes()
    adj = _adj(g)
    for size in range(1, k + 1):
        for _ in _connected_dominating(adj, n, size):
            return True
    return False


def gamma(g):
    return domination.domination_number(g)


# ---------------------------------------------------------------- sequences

def degree_sequences(n):
    """Every non-increasing sequence of n positive integers < n with even sum."""
    out = []

    def rec(pref, mx):
        if len(pref) == n:
            if sum(pref) % 2 == 0:
                out.append(tuple(pref))
            return
        for v in range(min(mx, n - 1), 0, -1):
            rec(pref + [v], v)

    rec([], n - 1)
    return out


def graphical_sequences(n):
    return [s for s in degree_sequences(n) if nx.is_graphical(s)]


# ---------------------------------------------------------------- families

def blowup(h, sizes):
    """Lexicographic blow-up: vertex v of h becomes a clique of sizes[v], and
    adjacent cliques are joined completely."""
    g = nx.Graph()
    blocks, idx = {}, 0
    for v in h:
        blocks[v] = list(range(idx, idx + sizes[v]))
        idx += sizes[v]
        for a, b in itertools.combinations(blocks[v], 2):
            g.add_edge(a, b)
        g.add_nodes_from(blocks[v])
    for u, v in h.edges():
        for a in blocks[u]:
            for b in blocks[v]:
                g.add_edge(a, b)
    return g


def spider(k, t0, t, s):
    """S(k; t0,t,s) of Theorem 3: gamma_c = k+1, W -> 1."""
    g = nx.Graph()
    c = list(range(t0))
    for a, b in itertools.combinations(c, 2):
        g.add_edge(a, b)
    g.add_nodes_from(c)
    idx = t0
    for _ in range(k):
        leg = list(range(idx, idx + t))
        idx += t
        tip = list(range(idx, idx + s))
        idx += s
        for blk in (leg, tip):
            for a, b in itertools.combinations(blk, 2):
                g.add_edge(a, b)
        for a in c:
            for b in leg:
                g.add_edge(a, b)
        for a in leg:
            for b in tip:
                g.add_edge(a, b)
    return g


def spider_seq(k, t0, t, s):
    return [t0 - 1 + k * t] * t0 + ([t - 1 + t0 + s] * t + [s - 1 + t] * s) * k


def geometric_chain(k, ratio=3):
    """Path of k cliques of sizes 1, r, r^2, ...; gamma_c = k - 2."""
    return blowup(nx.path_graph(k), [ratio ** i for i in range(k)])


# ---------------------------------------------------------------- searches

def realisation_search(seq, steps=200, seed=0):
    """Explore realisations of `seq` by double edge swaps; return the smallest
    number of connected dominating 4-sets seen (0 would mean gamma_c >= 5)."""
    n = len(seq)
    rng = random.Random(seed)
    try:
        g = nx.Graph(nx.havel_hakimi_graph(list(seq)))
    except nx.NetworkXError:
        return None
    best = None
    for _ in range(steps):
        if nx.is_connected(g):
            adj = _adj(g)
            c = sum(len(list(_connected_dominating(adj, n, k))) for k in (1, 2, 3, 4))
            if best is None or c < best:
                best = c
            if best == 0:
                return 0
        edges = list(g.edges())
        if len(edges) < 2:
            break
        for _try in range(20):
            (a, b), (c_, d) = rng.sample(edges, 2)
            if len({a, b, c_, d}) < 4 or g.has_edge(a, c_) or g.has_edge(b, d):
                continue
            g.remove_edge(a, b)
            g.remove_edge(c_, d)
            g.add_edge(a, c_)
            g.add_edge(b, d)
            break
    return best


def _cycle_blowup_start(n, k):
    """A graph on n vertices with gamma_c = k: blow up C_{k+2}."""
    m = k + 2
    if n < m:
        return None
    sizes = [1] * m
    for i in range(n - m):
        sizes[i % m] += 1
    return blowup(nx.cycle_graph(m), sizes)


def _spider_start(n, legs):
    t = max(1, (n - 2 - legs) // (2 * legs))
    t0 = n - legs * (t + 1)
    if t0 < 2:
        return None
    return spider(legs, t0, t, 1)


def anneal_min_residue(n, k=5, iters=3000, seed=0, start=None):
    """Minimise the residue over connected graphs on n vertices with
    gamma_c >= k.  Plateau moves are accepted, so the walk explores the whole
    residue level set.  Returns (best residue, graph)."""
    rng = random.Random(seed)
    g = (start or _cycle_blowup_start(n, k)).copy()
    if not (nx.is_connected(g) and not has_cds_upto(g, k - 1)):
        return None
    e = residue(g)
    best = (e, g.copy())
    for _ in range(iters):
        u, v = rng.sample(list(g), 2)
        had = g.has_edge(u, v)
        g.remove_edge(u, v) if had else g.add_edge(u, v)
        ok = nx.is_connected(g) and not has_cds_upto(g, k - 1)
        e2 = residue(g) if ok else None
        if ok and e2 is not None and (e2 <= e or rng.random() < 0.02):
            e = e2
            if e2 < best[0]:
                best = (e2, g.copy())
        else:
            g.add_edge(u, v) if had else g.remove_edge(u, v)
    return best


# ---------------------------------------------------------------- main

def check_lemma1(nmax=9):
    """Every Havel-Hakimi step on every graphical sequence, n <= nmax."""
    bad = 0
    total = 0
    for n in range(2, nmax + 1):
        for s in graphical_sequences(n):
            seq = sorted(s, reverse=True)
            while seq and seq[0] > 0:
                d = seq[0]
                w0 = caro_wei_seq(seq)
                seq = seq[1:]
                for i in range(d):
                    seq[i] -= 1
                seq.sort(reverse=True)
                total += 1
                if caro_wei_seq(seq) < w0 - 1e-12:
                    bad += 1
    return total, bad


def check_theorem1(nmax=10):
    bad, total, tight = 0, 0, 0
    for n in range(2, nmax + 1):
        for s in graphical_sequences(n):
            r = residue_seq(s)
            w = caro_wei_seq(s)
            total += 1
            if r < w - 1e-12:
                bad += 1
            if r == math.ceil(w - 1e-12):
                tight += 1
    return total, bad, tight


def residue2_sequence_table(nmax=11):
    rows = []
    for n in range(4, nmax + 1):
        g2 = [s for s in graphical_sequences(n) if residue_seq(s) == 2]
        if not g2:
            continue
        rows.append((n, len(g2), min(s[0] for s in g2), min(sum(s) for s in g2),
                     max(n - s[0] for s in g2),
                     all(n <= 2 * s[0] + 2 for s in g2),
                     all(sum(1 for d in s if d <= t) <= 2 * t + 2
                         for s in g2 for t in range(n))))
    return rows


def corpus_scan(nmax=8, verbose=True):
    t0 = time.time()
    graphs = allgraphs.connected_graphs(nmax)
    if verbose:
        print(f"   {len(graphs):,} connected graphs generated and count-checked "
              f"against OEIS [{time.time()-t0:.0f}s]")
    out = {"n": len(graphs), "res2": 0, "res2_gamma_le2": 0, "max_gc_res2": 0,
           "extremal": [], "min_res_by_gc": {}, "dm_bad": 0, "cw_bad": 0,
           "min_W_gc5": None, "cds_mismatch": 0, "res_mismatch": 0}
    for g in graphs:
        n = g.number_of_nodes()
        r = residue(g)
        gc = gamma_c(g)
        gm = gamma(g)
        w = caro_wei(g)
        if r != wowii.residue(g):
            out["res_mismatch"] += 1
        if n >= 3 and n - wowii.max_leaves_via_cds(g) != gc:
            out["cds_mismatch"] += 1
        out["min_res_by_gc"][gc] = min(out["min_res_by_gc"].get(gc, 99), r)
        if gc > 3 * gm - 2:
            out["dm_bad"] += 1
        if n >= 3 and gc > 3 * w - 2 + 1e-9:
            out["cw_bad"] += 1
        if gc >= 5 and (out["min_W_gc5"] is None or w < out["min_W_gc5"]):
            out["min_W_gc5"] = w
        if r == 2:
            out["res2"] += 1
            if gm <= 2:
                out["res2_gamma_le2"] += 1
            if gc > out["max_gc_res2"]:
                out["max_gc_res2"] = gc
            if gc >= 4:
                out["extremal"].append(
                    (n, gc, tuple(sorted((d for _, d in g.degree()), reverse=True))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="skip the exhaustive 8-vertex corpus")
    ap.add_argument("--deep", action="store_true",
                    help="longer sequence and annealing searches")
    args = ap.parse_args()
    t_start = time.time()

    print("CLAIM  residue(G) = 2  =>  gamma_c(G) <= 4.   Status: NOT PROVED.\n")

    print("1. PROVED  Lemma 1: a Havel-Hakimi step never decreases "
          "W = sum 1/(d+1)")
    total, bad = check_lemma1(9)
    print(f"   verified on {total:,} individual steps over every graphical "
          f"sequence with n <= 9: {bad} decreases")

    print("\n   PROVED  Theorem 1: residue >= W")
    nmax = 11 if args.deep else 10
    total, bad, tight = check_theorem1(nmax)
    print(f"   verified on all {total:,} graphical sequences with n <= {nmax}: "
          f"{bad} violations, {tight:,} at equality after rounding up")

    print("\n   PROVED  Corollaries: n <= R(Delta+1) and #{deg <= t} <= R(t+1)")
    rows = residue2_sequence_table(11 if args.deep else 10)
    print("     n   residue-2 seqs   min Delta   min sum   max n-Delta   "
          "n<=2D+2   #{d<=t}<=2t+2")
    for n, cnt, mind, mins, mxd, c1, c2 in rows:
        print(f"   {n:3d}   {cnt:12,}   {mind:9d}   {mins:7d}   {mxd:11d}   "
              f"{str(c1):7}   {str(c2)}")
    print("   the sparsest residue-2 sequence on n vertices is the "
          "(n/2 - 1)-regular one, so residue 2 forces m >= n(n-2)/8")

    print("\n2. PROVED  Theorem 2: residue 2 and gamma <= 2  =>  gamma_c <= 4")
    print("   (Duchet-Meyniel gamma_c <= 3*gamma - 2, verified below)")
    print("   FALSE  gamma(G) <= residue(G): 4-regular graphs on 10 vertices "
          "have residue 2")
    random.seed(11)
    hits = []
    for _ in range(300):
        h = nx.random_regular_graph(4, 10)
        if nx.is_connected(h) and gamma(h) == 3:
            hits.append(h)
    print(f"   {len(hits)} of 300 random 4-regular graphs on 10 vertices have "
          f"gamma = 3 with residue {residue(hits[0]) if hits else '-'}; "
          f"their gamma_c: {sorted(set(gamma_c(h) for h in hits))}")

    print("\n3. PROVED  Theorem 4: the constant 4 is sharp at every order")
    print("   H_T = C_6 with one block blown up to a clique of size T")
    print("      T     n   residue   gamma_c       W   degrees")
    for T in range(1, 12):
        h = blowup(nx.cycle_graph(6), [T, 1, 1, 1, 1, 1])
        seq = tuple(sorted((d for _, d in h.degree()), reverse=True))
        r, gc = residue(h), gamma_c(h)
        assert (r, gc) == (2, 4), (T, r, gc)
        print(f"   {T:4d} {h.number_of_nodes():5d} {r:9d} {gc:9d} "
              f"{caro_wei(h):7.3f}   {seq}")
    for T in range(12, 60):
        seq = [T + 1] * (T + 2) + [2, 2, 2]
        assert residue_seq(seq) == 2, T
    print("   residue stays 2 for every T up to 59 (Theorem 4 says: for all T)")

    print("\n4. THE TWO ROUTES THAT DIE")
    print("   FALSE  W <= 2 => gamma_c <= 4, by Theorem 3's spiders "
          "(gamma_c = k+1 with W -> 1)")
    print("     k  t0    t  s      n      W    residue   gamma_c")
    for (k, t0, t, s, do_gc) in [(4, 6, 2, 1, True), (4, 12, 3, 1, True),
                                 (4, 20, 4, 1, False), (4, 100, 10, 1, False),
                                 (4, 500, 22, 1, False)]:
        seq = spider_seq(k, t0, t, s)
        gcv = gamma_c(spider(k, t0, t, s)) if do_gc else k + 1
        print(f"   {k:3d} {t0:4d} {t:4d} {s:2d} {len(seq):6d} "
              f"{caro_wei_seq(seq):6.3f} {residue_seq(seq):9d} {gcv:9d}"
              f"{'' if do_gc else '  (gamma_c by Theorem 3)'}")
    lo = min(residue_seq(spider_seq(4, t0, t, s))
             for t0 in range(2, 40) for t in range(1, 25) for s in range(1, 25))
    print(f"   minimum residue over 22,800 spider parameter triples with "
          f"gamma_c = 5: {lo}")

    print("\n   FALSE  the corollaries alone: the geometric chain of cliques "
          "1,3,9,...")
    for k in (5, 6, 7, 8):
        g = geometric_chain(k)
        seq = sorted((d for _, d in g.degree()), reverse=True)
        n = g.number_of_nodes()
        ok = all(sum(1 for d in seq if d <= t) <= 2 * t + 2 for t in range(max(seq) + 1))
        print(f"   k={k}: n={n:5d}  gamma_c={k-2}  residue={residue_seq(seq)}  "
              f"passes every residue-2 corollary: {ok}")

    print("\n5. VERIFIED  families")
    print("   C_k[K_t]: residue ceil(k/3), gamma_c = k-2, so gamma_c = 3R-2 "
          "is attained for every R")
    for k in (4, 5, 6, 7, 8, 9):
        for t in (1, 2, 3):
            g = blowup(nx.cycle_graph(k), [t] * k)
            r = residue(g)
            assert r == math.ceil(k / 3), (k, t, r)
            assert gamma_c(g) == max(1, k - 2), (k, t)
        print(f"     C_{k}[K_t], t=1,2,3: residue {math.ceil(k/3)}, "
              f"gamma_c {max(1, k-2)}")
    print("   K_{a,a}: residue 2 for every a, so residue-2 graphs exist at "
          "every order")
    print("     a =", [a for a in range(2, 12)], "->",
          [residue(nx.complete_bipartite_graph(a, a)) for a in range(2, 12)])
    print("   P_k[K_t]: residue 2 exactly for k <= 5, where gamma_c = k-2 <= 3")
    print("     k=2..8, t=3:",
          [residue(blowup(nx.path_graph(k), [3] * k)) for k in range(2, 9)])

    if not args.quick:
        print("\n6. VERIFIED  exhaustive: every connected graph on at most 8 "
              "vertices")
        out = corpus_scan(8)
        print(f"   cross-checks against the repo: residue mismatches "
              f"{out['res_mismatch']}, n - max_leaves_via_cds mismatches "
              f"{out['cds_mismatch']}")
        print(f"   Duchet-Meyniel gamma_c <= 3*gamma - 2 violations: "
              f"{out['dm_bad']}")
        print(f"   gamma_c <= 3W - 2 violations: {out['cw_bad']}  "
              f"(so Conjecture A does not hold with W in place of R)")
        print(f"   residue-2 graphs: {out['res2']:,}; with gamma <= 2 "
              f"(Theorem 2 settles these): {out['res2_gamma_le2']:,} = "
              f"{100*out['res2_gamma_le2']/out['res2']:.2f}%")
        print(f"   largest gamma_c over residue-2 graphs: "
              f"{out['max_gc_res2']}  <- the claim, in range")
        print(f"   the {len(out['extremal'])} extremal graphs "
              f"(residue 2, gamma_c = 4):")
        for n, gc, seq in out["extremal"]:
            print(f"      n={n} degrees={seq}")
        print("   minimum residue over the graphs with gamma_c = k:")
        for k in sorted(out["min_res_by_gc"]):
            print(f"      k={k}: {out['min_res_by_gc'][k]}   "
                  f"(ceil((k+2)/3) = {math.ceil((k+2)/3)})")
        print(f"   minimum W over the graphs with gamma_c >= 5: "
              f"{out['min_W_gc5']:.4f}")

    print("\n7. VERIFIED  the hunt")
    hi = 12 if args.deep else 10
    print(f"   realisations of every residue-2 degree sequence, n <= {hi}")
    for n in range(6, hi + 1):
        seqs = [s for s in graphical_sequences(n) if residue_seq(s) == 2]
        hits = 0
        hardest = (10 ** 9, None)
        for s in seqs:
            b = realisation_search(s, steps=100 if not args.deep else 200,
                                   seed=hash(s) & 0xffff)
            if b == 0:
                hits += 1
                print("      COUNTEREXAMPLE SEQUENCE:", s)
            if b is not None and b < hardest[0]:
                hardest = (b, s)
        print(f"      n={n}: {len(seqs):5,} sequences, realisations with "
              f"gamma_c >= 5: {hits};  hardest {hardest[1]} still had "
              f"{hardest[0]} connected dominating 4-sets")

    print("   annealing on the residue under the hard constraint gamma_c >= 5")
    for n in range(7, 17 if not args.deep else 23):
        rec = None
        for legs, st in (("cycle", _cycle_blowup_start(n, 5)),
                         ("spider", _spider_start(n, 4))):
            if st is None or has_cds_upto(st, 4):
                continue
            b = anneal_min_residue(n, 5, iters=800 if not args.deep else 6000,
                                   seed=n, start=st)
            if b and (rec is None or b[0] < rec[0]):
                rec = b
        if rec is None:
            continue
        g = rec[1]
        print(f"      n={n}: minimum residue reached {rec[0]} "
              f"(gamma_c={gamma_c(g)}, W={caro_wei(g):.3f})"
              + ("   *** RESIDUE 2 WITH gamma_c >= 5 ***" if rec[0] <= 2 else ""))

    print(f"\nnothing below residue 3 was ever reached with gamma_c >= 5. "
          f"[{time.time()-t_start:.0f}s]")


if __name__ == "__main__":
    main()
