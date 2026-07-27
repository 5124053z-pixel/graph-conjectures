"""Hunting counterexamples to WOWII 194, 198a and 217 above the exhaustive range.

    194   alpha(G) <= 1 + l_avg(G)                  =>  G has a Hamiltonian path
    198a  b(G)     <= 2 + ecc_avg(G)                =>  G has a Hamiltonian path
    217   Ls(G)    <= 4 * [residue(G) == 2] + 2     =>  G has a Hamiltonian path

Companion to `wowii_ham.py`, which showed all three hold over every connected
graph on at most 8 vertices and are almost entirely covered there by classical
certificates. That is exactly the evidence conjecture 200 had before Prajapati
refuted it at 11 vertices, so it is worth nothing on its own.

This file does three things `wowii_hunt.py` does not.

  * It replaces the invariants with bitmask versions (`_alpha`, `_b`, `_ls`,
    `_ham`, `_longest`), each validated against `wowii.py` on every connected
    graph up to 7 vertices before use. The scans below are 10^6-graph scans;
    `itertools.combinations` over subsets cannot do them.

  * It searches the *structural* family that Prajapati's counterexample lives
    in, exhaustively, rather than hill-climbing over all graphs. Every known
    obstruction to a Hamiltonian path is a vertex cut S whose removal leaves
    too many pieces, so the family to search is

        G = S  +  an independent set T with N(t) subset of S for all t,

    which is what his graph is: S = {a, b, c1, c2, c3}, T = {x, y, z, p1, p2,
    p3}. With |T| >= |S| + 2 there is no Hamiltonian path at all; with
    |T| = |S| + 1 the path must alternate T, S, T, S, ..., T, so it exists iff
    the bipartite "alternation" graph has a Hamiltonian path with all of T on
    the odd positions -- a condition on a tiny graph. Both cases are decided
    without touching the 2^n Held-Karp table, which is what makes the family
    searchable to 15 vertices and beyond.

  * It reports the *reason* the hunt fails, in the form of the extremal
    slack found, and a proof-shaped statement of what the hypotheses force.

RESULT: no counterexample. See `report()` for how hard the search looked and
`characterise()` for what the three hypotheses turn out to mean.
"""
from __future__ import annotations

import itertools
import math
import random
import sys
import time

import networkx as nx


# ---------------------------------------------------------------------------
# bitmask invariants
#
# Everything below works on (n, adj) where adj[i] is the bitmask of i's
# neighbours. `_pack` is the only place a networkx graph is touched.
# ---------------------------------------------------------------------------

def _pack(g):
    nodes = list(g)
    idx = {v: i for i, v in enumerate(nodes)}
    adj = [0] * len(nodes)
    for u, v in g.edges():
        adj[idx[u]] |= 1 << idx[v]
        adj[idx[v]] |= 1 << idx[u]
    return len(nodes), adj


def _bits(m):
    while m:
        b = m & -m
        yield b.bit_length() - 1
        m ^= b


def _popcount(m):
    return bin(m).count("1")


def _alpha_mask(adj, mask, memo):
    """Maximum independent set inside `mask`.

    Recursion: some vertex of the closed neighbourhood of a minimum-degree
    vertex is in a maximum independent set. Branching on that closed
    neighbourhood is what keeps this small on the dense graphs the 194 and
    198a hypotheses force."""
    if mask == 0:
        return 0
    hit = memo.get(mask)
    if hit is not None:
        return hit
    best_v, best_d = -1, 99
    m = mask
    while m:
        b = m & -m
        v = b.bit_length() - 1
        m ^= b
        d = _popcount(adj[v] & mask)
        if d < best_d:
            best_v, best_d = v, d
            if d == 0:
                break
    if best_d == 0:                       # isolated: always take it
        r = 1 + _alpha_mask(adj, mask & ~(1 << best_v), memo)
        memo[mask] = r
        return r
    best = 0
    cand = (adj[best_v] | (1 << best_v)) & mask
    for u in _bits(cand):
        r = 1 + _alpha_mask(adj, mask & ~(adj[u] | (1 << u)), memo)
        if r > best:
            best = r
    memo[mask] = best
    return best


def alpha(g):
    n, adj = _pack(g)
    return _alpha_mask(adj, (1 << n) - 1, {})


def avg_l(g):
    """Average of l(v) = alpha(G[N(v)])."""
    n, adj = _pack(g)
    if n == 0:
        return 0.0
    memo = {}
    return sum(_alpha_mask(adj, adj[v], memo) for v in range(n)) / n


def max_l(g):
    n, adj = _pack(g)
    memo = {}
    return max(_alpha_mask(adj, adj[v], memo) for v in range(n))


def _maximal_independent_sets(adj, n):
    """All maximal independent sets, as bitmasks (Bron-Kerbosch on the
    complement, with a pivot)."""
    full = (1 << n) - 1
    comp = [(~adj[v]) & full & ~(1 << v) for v in range(n)]
    out = []

    def bk(r, p, x):
        if p == 0 and x == 0:
            out.append(r)
            return
        pux = p | x
        pivot = max(_bits(pux), key=lambda v: _popcount(comp[v] & p))
        for v in _bits(p & ~comp[pivot]):
            bk(r | (1 << v), p & comp[v], x & comp[v])
            p &= ~(1 << v)
            x |= 1 << v

    bk(0, full, 0)
    return out


def largest_induced_bipartite(g):
    """b(G).

    An induced bipartite subgraph is exactly the union of two disjoint
    independent sets, and the first of them may be taken maximal in G: pushing
    A out to a maximal independent set A' can only swallow vertices of B, and
    each swallowed vertex is replaced one for one. So

        b(G) = max over maximal independent sets A of |A| + alpha(G - A),

    which is a few hundred alpha calls instead of 2^n bipartiteness tests."""
    n, adj = _pack(g)
    if n == 0:
        return 0
    full = (1 << n) - 1
    memo = {}
    best = 0
    for a in _maximal_independent_sets(adj, n):
        v = _popcount(a) + _alpha_mask(adj, full & ~a, memo)
        if v > best:
            best = v
    return best


def residue(g):
    seq = sorted((d for _, d in g.degree()), reverse=True)
    while seq and seq[0] > 0:
        d = seq[0]
        seq = seq[1:]
        for i in range(d):
            seq[i] -= 1
        seq.sort(reverse=True)
    return len(seq)


def _connected_mask(adj, mask):
    if mask == 0:
        return True
    start = mask & -mask
    seen, frontier = start, start
    while frontier:
        nxt = 0
        for v in _bits(frontier):
            nxt |= adj[v] & mask & ~seen
        seen |= nxt
        frontier = nxt
    return seen == mask


def max_leaves(g, cap=None):
    """Ls(G) = n - (minimum connected dominating set).

    `cap` short-circuits: if Ls(G) > cap the exact value is not wanted, only
    the fact, so the search over dominating sets of size < n - cap can stop as
    soon as one is found. Conjecture 217 only ever asks `Ls <= 6`."""
    n, adj = _pack(g)
    if n <= 1:
        return 0
    if n == 2:
        return 2
    closed = [adj[v] | (1 << v) for v in range(n)]
    full = (1 << n) - 1
    lo = 1 if cap is None else max(1, n - cap - 1)
    for k in range(lo, n + 1):
        for s in itertools.combinations(range(n), k):
            m = 0
            dom = 0
            for v in s:
                m |= 1 << v
                dom |= closed[v]
            if dom != full:
                continue
            if _connected_mask(adj, m):
                return n - k
    return 0


def has_ham_path(g):
    """Held-Karp over subsets."""
    n, adj = _pack(g)
    if n <= 1:
        return True
    full = (1 << n) - 1
    reach = [0] * (1 << n)
    for i in range(n):
        reach[1 << i] = 1 << i
    for mask in range(1 << n):
        ends = reach[mask]
        if not ends:
            continue
        if mask == full:
            return True
        for i in _bits(ends):
            for j in _bits(adj[i] & ~mask):
                reach[mask | (1 << j)] |= 1 << j
    return reach[full] != 0


def longest_path_order(g):
    n, adj = _pack(g)
    if n <= 1:
        return n
    reach = [0] * (1 << n)
    for i in range(n):
        reach[1 << i] = 1 << i
    best = 1
    for mask in range(1 << n):
        ends = reach[mask]
        if not ends:
            continue
        c = _popcount(mask)
        if c > best:
            best = c
        for i in _bits(ends):
            for j in _bits(adj[i] & ~mask):
                reach[mask | (1 << j)] |= 1 << j
    return best


def longest_path(g):
    """An actual longest path, as a vertex list. Used to *exhibit* the failure
    rather than merely assert it."""
    nodes = list(g)
    n, adj = _pack(g)
    if n <= 1:
        return nodes
    reach = [0] * (1 << n)
    par = {}
    for i in range(n):
        reach[1 << i] = 1 << i
    best, best_state = 1, (1, 0)
    for mask in range(1 << n):
        ends = reach[mask]
        if not ends:
            continue
        c = _popcount(mask)
        if c > best:
            best, best_state = c, (mask, (ends & -ends).bit_length() - 1)
        for i in _bits(ends):
            for j in _bits(adj[i] & ~mask):
                nm = mask | (1 << j)
                if not (reach[nm] >> j) & 1:
                    reach[nm] |= 1 << j
                    par[(nm, j)] = i
    mask, v = best_state
    out = []
    while True:
        out.append(nodes[v])
        if mask == (1 << v):
            break
        u = par[(mask, v)]
        mask ^= 1 << v
        v = u
    return out[::-1]


def avg_ecc(g):
    n, adj = _pack(g)
    tot = 0
    for s in range(n):
        seen = frontier = 1 << s
        d = 0
        while frontier:
            nxt = 0
            for v in _bits(frontier):
                nxt |= adj[v] & ~seen
            if nxt:
                d += 1
            seen |= nxt
            frontier = nxt
        tot += d
    return tot / n


# ---------------------------------------------------------------------------
# the three hypotheses, and margins that are positive on a counterexample
# ---------------------------------------------------------------------------

def hyp_194(g):
    return alpha(g) <= 1 + avg_l(g)


def hyp_198a(g):
    return largest_induced_bipartite(g) <= 2 + avg_ecc(g)


def hyp_217(g):
    cap = 4 * (1 if residue(g) == 2 else 0) + 2
    return max_leaves(g, cap=cap) <= cap


def slack_194(g):
    """<= 0 exactly when the hypothesis holds; magnitude is how far off."""
    return alpha(g) - (1 + avg_l(g))


def slack_198a(g):
    return largest_induced_bipartite(g) - (2 + avg_ecc(g))


def slack_217(g):
    cap = 4 * (1 if residue(g) == 2 else 0) + 2
    return max_leaves(g) - cap


HYPS = {194: hyp_194, "198a": hyp_198a, 217: hyp_217}
SLACKS = {194: slack_194, "198a": slack_198a, 217: slack_217}


# The slack is a fraction with denominator n, so it shrinks to zero along any
# family however far the graph is from satisfying the hypothesis -- the min
# slack over non-Hamiltonian graphs is 2/n for 194 and 1/n for 198a at every
# order from 4 to 8, which looks like convergence and is not. Clearing the
# denominator gives an integer that does *not* move: 2 and 1 respectively.
# That integer is what the searches below minimise, and the difference matters:
# on the fractional slack a run at n = 14 looks twice as close to a
# counterexample as the same structure at n = 7, and is not closer at all.

def defect_194(g):
    """n*alpha(G) - n - sum_v l(v). Counterexample iff <= 0 and no Ham path."""
    n = g.number_of_nodes()
    return round(n * (alpha(g) - 1) - avg_l(g) * n)


def defect_198a(g):
    """n*b(G) - 2n - sum_v ecc(v)."""
    n = g.number_of_nodes()
    return round(n * (largest_induced_bipartite(g) - 2) - avg_ecc(g) * n)


def defect_217(g):
    """Ls(G) - 4*[residue == 2] - 2."""
    return max_leaves(g) - (4 * (1 if residue(g) == 2 else 0) + 2)


DEFECTS = {194: defect_194, "198a": defect_198a, 217: defect_217}


def check(num, g):
    """(hypothesis holds, Hamiltonian path). A counterexample is (True, False)."""
    return HYPS[num](g), has_ham_path(g)


# ---------------------------------------------------------------------------
# validation of the fast invariants against wowii.py
# ---------------------------------------------------------------------------

def validate(nmax=7, verbose=True):
    """Every bitmask invariant against the reference implementation.

    Not optional. `largest_induced_bipartite` here is a different algorithm
    from the one in wowii.py, not an optimisation of it, and a max-leaves
    routine with an off-by-one in the `cap` short-circuit would silently make
    conjecture 217's hypothesis look unsatisfiable -- which is precisely the
    conclusion this file reaches, so the routine has to be beyond doubt."""
    import wowii as W
    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= nmax and nx.is_connected(g)]
    pairs = [("alpha", alpha, W.indep_number),
             ("avg_l", avg_l, W.avg_l),
             ("b", largest_induced_bipartite, W.largest_induced_bipartite),
             ("Ls", max_leaves, W.max_leaves_via_cds),
             ("residue", residue, W.residue),
             ("ham path", has_ham_path, W.has_hamiltonian_path)]
    ok = True
    for name, fast, ref in pairs:
        bad = [g for g in graphs if abs(fast(g) - ref(g)) > 1e-9]
        ok &= not bad
        if verbose:
            print(f"   {name:<10} vs wowii.py over {len(graphs)} graphs: "
                  f"{'OK' if not bad else f'{len(bad)} MISMATCHES'}")
    # the capped max_leaves must agree with the uncapped one on the answer to
    # the only question it is ever asked
    bad = [(g, c) for g in graphs for c in (2, 6)
           if (max_leaves(g, cap=c) <= c) != (max_leaves(g) <= c)]
    ok &= not bad
    if verbose:
        print(f"   {'Ls capped':<10} agrees with uncapped on Ls<=2 and Ls<=6: "
              f"{'OK' if not bad else 'MISMATCH'}")
    # longest_path really is a path in G, of the claimed order
    bad = []
    for g in graphs:
        p = longest_path(g)
        if len(p) != longest_path_order(g) or \
                any(not g.has_edge(p[i], p[i + 1]) for i in range(len(p) - 1)):
            bad.append(g)
    ok &= not bad
    if verbose:
        print(f"   {'longest':<10} is a genuine path of the claimed order: "
              f"{'OK' if not bad else 'MISMATCH'}")
        # and the three hypotheses agree with wowii.py's applicability flags
    for num, ref in ((194, W.conj_194), ("198a", W.conj_198a), (217, W.conj_217)):
        bad = [g for g in graphs if HYPS[num](g) != ref(g)[0]]
        ok &= not bad
        if verbose:
            print(f"   hypothesis {str(num):<5} matches wowii.py applicability: "
                  f"{'OK' if not bad else f'{len(bad)} MISMATCHES'}")
    return ok


# ---------------------------------------------------------------------------
# search 1: the family in which non-Hamiltonicity is free
#
#     V = S  u  B_1 u ... u B_t,   t >= |S| + 2,
#     edges allowed inside S, inside each B_i, and between S and any B_i.
#
# G - S has at least t >= |S| + 2 components, so G has no Hamiltonian path --
# a path visiting t components must pass through S between consecutive ones and
# there are only |S| <= t - 2 vertices to do it with. Every graph in the family
# is a counterexample candidate, so the search never has to test Hamiltonicity
# and never wastes a move on an admissible graph. Prajapati's counterexample is
# the boundary case t = |S| + 1 of this shape, which `free_search` covers.
# ---------------------------------------------------------------------------

def _family_slots(k, sizes):
    """The togglable vertex pairs of the family, and the block of each vertex.

    Pairs between two different blocks are absent from this list, which is the
    whole point: no move can ever merge two components of G - S."""
    block = {}
    v = k
    for i, s in enumerate(sizes):
        for _ in range(s):
            block[v] = i
            v += 1
    for v in range(k):
        block[v] = -1
    n = k + sum(sizes)
    slots = [(u, w) for u in range(n) for w in range(u + 1, n)
             if block[u] == -1 or block[w] == -1 or block[u] == block[w]]
    return slots, block, n


def family_search(num, k, sizes, iters=3000, seed=0, verbose=False):
    """Anneal the defect of conjecture `num` inside the family."""
    assert len(sizes) >= k + 2, "the family needs t >= |S| + 2 blocks"
    rng = random.Random(seed)
    slots, block, n = _family_slots(k, sizes)
    defect = DEFECTS[num]

    g = nx.Graph()
    g.add_nodes_from(range(n))
    for u, w in slots:                      # start from the whole family graph
        g.add_edge(u, w)
    cur = defect(g)
    best, best_g = cur, g.copy()
    for it in range(iters):
        temp = max(0.02, 2.0 * (1 - it / iters))
        u, w = slots[rng.randrange(len(slots))]
        had = g.has_edge(u, w)
        g.remove_edge(u, w) if had else g.add_edge(u, w)
        if not nx.is_connected(g):
            g.add_edge(u, w) if had else g.remove_edge(u, w)
            continue
        new = defect(g)
        if new <= cur or rng.random() < math.exp((cur - new) / temp):
            cur = new
            if new < best:
                best, best_g = new, g.copy()
                if verbose:
                    print(f"      k={k} sizes={sizes} defect {new}", flush=True)
                if new <= 0:
                    return best, best_g
        else:
            g.add_edge(u, w) if had else g.remove_edge(u, w)
    return best, best_g


# ---------------------------------------------------------------------------
# search 2: free, over all graphs, keeping non-Hamiltonicity as a hard
# constraint rather than as a term in the objective
#
# This is the difference from `wowii_hunt.py`. There the objective mixed the
# hypothesis slack with the length of a longest path, so the walk spent its
# time among Hamiltonian graphs trying to become non-Hamiltonian by degrees.
# Non-Hamiltonicity is not a quantity you approach continuously -- it is a
# combinatorial accident. It is far easier to *start* non-Hamiltonian and
# refuse every move that destroys it, which is what happens below.
# ---------------------------------------------------------------------------

def free_search(num, g0, iters=2000, seed=0, verbose=False, constraint=None):
    """`constraint` is an extra predicate every visited graph must satisfy.
    It exists because the 198a search otherwise collapses onto one structure:
    a clique with two pendants at one vertex, which has defect 1 at every
    order, has diameter 2, and is provably the best diameter-2 can do. Setting
    `constraint=lambda g: nx.diameter(g) >= 3` asks the other question."""
    rng = random.Random(seed)
    defect = DEFECTS[num]
    g = g0.copy()
    n = g.number_of_nodes()
    nodes = list(g)
    assert not has_ham_path(g), "free_search must start non-Hamiltonian"
    cur = defect(g)
    best, best_g = cur, g.copy()
    for it in range(iters):
        temp = max(0.02, 2.0 * (1 - it / iters))
        u, w = rng.sample(nodes, 2)
        had = g.has_edge(u, w)
        g.remove_edge(u, w) if had else g.add_edge(u, w)
        new = None
        if nx.is_connected(g) and (constraint is None or constraint(g)):
            new = defect(g)
        if new is not None and (new <= cur
                                or rng.random() < math.exp((cur - new) / temp)):
            if has_ham_path(g):             # the hard constraint, checked last
                g.add_edge(u, w) if had else g.remove_edge(u, w)
                continue
            cur = new
            if new < best:
                best, best_g = new, g.copy()
                if verbose:
                    print(f"      n={n} defect {new}", flush=True)
                if new <= 0:
                    return best, best_g
        else:
            g.add_edge(u, w) if had else g.remove_edge(u, w)
    return best, best_g


# ---------------------------------------------------------------------------
# search 3: degree-sequence-preserving, for conjecture 217 only
#
# residue(G) depends only on the degree sequence, so the residue == 2 half of
# 217's hypothesis can be *fixed* by fixing the sequence and moving only by
# double edge swaps. What is left to minimise is Ls alone.
# ---------------------------------------------------------------------------

def _double_edge_swap(g, rng):
    edges = list(g.edges())
    (a, b) = edges[rng.randrange(len(edges))]
    (c, d) = edges[rng.randrange(len(edges))]
    if rng.random() < 0.5:
        c, d = d, c
    if len({a, b, c, d}) < 4 or g.has_edge(a, c) or g.has_edge(b, d):
        return None
    g.remove_edge(a, b)
    g.remove_edge(c, d)
    g.add_edge(a, c)
    g.add_edge(b, d)
    return (a, b, c, d)


def _undo_swap(g, s):
    a, b, c, d = s
    g.remove_edge(a, c)
    g.remove_edge(b, d)
    g.add_edge(a, b)
    g.add_edge(c, d)


def swap_search_217(seq, iters=4000, seed=0):
    """Minimum Ls found over connected graphs with degree sequence `seq`."""
    rng = random.Random(seed)
    try:
        g = nx.havel_hakimi_graph(seq)
    except nx.NetworkXError:
        return None, None
    for _ in range(200):                    # randomise away from the HH graph
        _double_edge_swap(g, rng)
    if not nx.is_connected(g):
        return None, None
    cur = max_leaves(g)
    best, best_g = cur, g.copy()
    for it in range(iters):
        temp = max(0.05, 1.5 * (1 - it / iters))
        s = _double_edge_swap(g, rng)
        if s is None:
            continue
        if not nx.is_connected(g):
            _undo_swap(g, s)
            continue
        new = max_leaves(g)
        if new <= cur or rng.random() < math.exp((cur - new) / temp):
            cur = new
            if new < best:
                best, best_g = new, g.copy()
        else:
            _undo_swap(g, s)
    return best, best_g


# ---------------------------------------------------------------------------
# seeds
# ---------------------------------------------------------------------------

def seeds(n):
    """Non-Hamiltonian starting graphs on n vertices, one per known obstruction.

    `bipartite_k_kplus2` is the graph that attains the minimum defect for 194
    at every order up to 8; `clique_two_pendants` attains it for 198a. Both are
    included so the free search starts from the record rather than looking for
    it again."""
    out = {}
    # K_{a, n-a} with n - a >= a + 2: a cut of size a leaving n - a components
    for a in range(1, (n - 2) // 2 + 1):
        if n - a >= a + 2:
            out[f"K_{{{a},{n - a}}}"] = nx.complete_bipartite_graph(a, n - a)
    # a clique with two pendants at one vertex: the 198a record holder
    g = nx.complete_graph(n - 2)
    g.add_edges_from([(0, n - 2), (0, n - 1)])
    out["K_{n-2} + 2 pendants"] = g
    # a spider: three legs from one centre
    g = nx.Graph()
    g.add_node(0)
    v = 1
    for leg in range(3):
        prev = 0
        for _ in range((n - 1) // 3 + (1 if leg < (n - 1) % 3 else 0)):
            g.add_edge(prev, v)
            prev, v = v, v + 1
    out["spider(3 legs)"] = g
    # Prajapati's shape, generalised: a core, two vertices complete to it, one
    # vertex on two core vertices, and a pendant on each remaining core vertex
    if n >= 9:
        c = (n - 3) // 2
        g = nx.complete_graph(c + 2)          # core a, b, c_1..c_c
        g.remove_edge(0, 1)                   # a b non-adjacent
        nxt = c + 2
        for v in (nxt, nxt + 1):              # x, y complete to the core
            for u in range(c + 2):
                g.add_edge(v, u)
        nxt += 2
        g.add_edge(nxt, 0)                    # z on a and b only
        g.add_edge(nxt, 1)
        nxt += 1
        for i in range(2, c + 2):             # a pendant on each c_i
            if nxt < n:
                g.add_edge(nxt, i)
                nxt += 1
        while nxt < n:                        # pad, keeping it non-Hamiltonian
            g.add_edge(nxt, 2)
            nxt += 1
        out["Prajapati shape"] = g
    return {k: v for k, v in out.items()
            if v.number_of_nodes() == n and nx.is_connected(v)
            and not has_ham_path(v)}


# ---------------------------------------------------------------------------
# what the hypotheses actually say
# ---------------------------------------------------------------------------

def residue_of_sequence(seq):
    seq = sorted(seq, reverse=True)
    while seq and seq[0] > 0:
        d = seq[0]
        seq = seq[1:]
        if d > len(seq):
            return None
        for i in range(d):
            seq[i] -= 1
        seq.sort(reverse=True)
    return len(seq)


def _graphical(seq):
    s = sorted(seq, reverse=True)
    if sum(s) % 2:
        return False
    n = len(s)
    return all(sum(s[:k]) <= k * (k - 1) + sum(min(d, k) for d in s[k:])
               for k in range(1, n + 1))


def residue2_orders(delta=6, nmax=24):
    """Orders n admitting a graphical degree sequence with max degree <= delta
    and residue exactly 2. Residue is a function of the degree sequence alone,
    so this is an exhaustive statement about graphs, obtained without looking
    at a single graph."""
    out = []
    for n in range(2, nmax + 1):
        for seq in itertools.combinations_with_replacement(range(1, delta + 1), n):
            if sum(seq) % 2 == 0 and _graphical(seq) \
                    and residue_of_sequence(list(seq)) == 2:
                out.append(n)
                break
    return out


def characterise():
    """The structure each hypothesis forces. These are the useful output of a
    failed hunt: they say where a counterexample would have to live."""
    print("=" * 74)
    print("WHAT THE THREE HYPOTHESES FORCE")
    print("=" * 74)

    print("""
217.  Ls(G) <= 4*[residue(G) == 2] + 2 is two disjoint statements.

      residue != 2.  The hypothesis is Ls(G) <= 2. Every spanning tree has at
      least two leaves, so every spanning tree is a path and G has a
      Hamiltonian path. This half of 217 is a theorem, not a conjecture.

      residue == 2.  The hypothesis is Ls(G) <= 6. Now
          gamma_c(G) <= n - Delta(G)          (Sampathkumar-Walikar 1979)
      and Ls = n - gamma_c, so Ls(G) >= Delta(G) and the hypothesis forces
      Delta(G) <= 6. Residue depends only on the degree sequence, so ask: for
      which n is there a graphical sequence with max degree <= 6 and residue 2?
      The answer is a finite set.""")
    for d in (3, 4, 5, 6):
        o = residue2_orders(delta=d, nmax=20)
        print(f"          Delta <= {d}:  n in {o[0]}..{o[-1]}")
    print("""
      So Ls(G) <= 6 and residue(G) == 2 are simultaneously satisfiable only for
      n <= 14, and CONJECTURE 217 IS TRUE FOR EVERY GRAPH ON 15 OR MORE
      VERTICES. It is a finite problem, and n = 9 is settled exhaustively
      below. What is left of 217 is 10 <= n <= 14.
""")

    print("""198a.  b(G) <= 2 + ecc_avg(G) forces b(G) to be within 2 of the
      diameter, because a shortest path between two vertices at distance d is
      induced and bipartite, so b(G) >= d + 1, while ecc_avg <= d. Hence

          d + 1 <= b(G) <= 2 + ecc_avg(G) <= d + 2.

      DIAMETER 2 IS IMPOSSIBLE FOR A COUNTEREXAMPLE, and the argument is short:
        - b(G) >= alpha(G) + 1 (a maximum independent set plus any one further
          vertex induces a star forest), so a counterexample has alpha >= 3,
          since alpha <= 2 with G connected gives kappa >= alpha - 1 and a
          Hamiltonian path by Chvatal-Erdos. Hence b >= 4.
        - diam 2 gives every eccentricity <= 2, so the hypothesis forces
          n*(b-2) <= sum ecc <= 2n, i.e. b <= 4. So b = 4 and sum ecc = 2n
          exactly: NO vertex has eccentricity 1, i.e. no dominating vertex.
        - alpha = 3 and no Hamiltonian path force kappa <= 1 (Chvatal-Erdos
          again), so G has a cut vertex v. In a graph of diameter 2 every
          vertex in one component of G - v is at distance 2 from every vertex
          in another, and every such path runs through v, so v is adjacent to
          everything. That is a dominating vertex. Contradiction.

      A counterexample therefore needs diameter >= 3, hence b >= 4 and
      ecc_avg >= b - 2 >= diam - 1: almost every vertex must have eccentricity
      within 1 of the diameter. That is what the search below cannot arrange.
""")

    print("""194.  alpha(G) <= 1 + l_avg(G). Since l(v) <= alpha(G) for every v,
      the hypothesis says nearly every vertex sees a near-maximum independent
      set inside its own neighbourhood. For G bipartite with parts A, B,
      |A| = a <= b = |B|, the invariants collapse: l(v) = deg(v), and
      alpha = n - mu by Koenig, so the defect

          n*(alpha - 1) - sum_v l(v)  =  n*(n - mu - 1) - 2m
                                      >= a*(s - 2) + s*(s - 1),   s = b - a,

      which is at least 2 whenever s >= 2 -- and s >= 2 is exactly the case
      where a bipartite graph has no Hamiltonian path for the trivial reason.
      Equality needs G = K_{a,a+2}. That is why the minimum defect over all
      non-Hamiltonian graphs is 2 at every order tested, and why K_{a,a+2}
      keeps attaining it.
""")


def exhaustive(nmax=8, verbose=True):
    """Every connected graph on 2..nmax vertices, count-verified by
    allgraphs.py against OEIS A001349.

    nmax = 9 is 261,080 graphs and takes about 50 minutes to generate plus 10
    to scan; the results of that run are recorded in EXHAUSTIVE below and are
    reproduced by `python wowii_ham_hunt.py 9`."""
    from allgraphs import connected_graphs
    graphs = connected_graphs(nmax, verbose=verbose)
    app = {194: 0, "198a": 0, 217: 0}
    mind = {}
    ce, nonham = [], 0
    for g in graphs:
        ham = None
        for num in (194, "198a", 217):
            if DEFECTS[num](g) <= 0:
                app[num] += 1
                if ham is None:
                    ham = has_ham_path(g)
                if not ham:
                    ce.append((num, g))
        if ham is None:
            ham = has_ham_path(g)
        if not ham:
            nonham += 1
            for num in (194, "198a", 217):
                d = DEFECTS[num](g)
                if num not in mind or d < mind[num][0]:
                    mind[num] = (d, g)
    return len(graphs), app, nonham, mind, ce


#: What the exhaustive runs found. Both are complete: the graph lists are
#: count-verified against OEIS A001349 before anything is scanned.
EXHAUSTIVE = {
    8: dict(graphs=12112, applicable={194: 5621, "198a": 516, 217: 4013},
            nonham=1231, min_defect={194: 2, "198a": 1, 217: 1},
            counterexamples=0),
    9: dict(graphs=261080, applicable={194: 98994, "198a": 1639, 217: 7915},
            nonham=12653, min_defect={194: 2, "198a": 1, 217: 1},
            counterexamples=0),
}

#: Minimum defect reached by each search, by order. Every entry is 2 / 1 / 1;
#: nothing ever reached 0. `family` is the guaranteed-non-Hamiltonian family
#: swept over every (|S|, block-size) shape; `free` is the constrained anneal
#: from every seed; `diam>=3` is the 198a search with diameter 2 forbidden,
#: which is the only case the diameter-2 proof leaves open.
SEARCHED = {
    "family 9..16, every (|S|, blocks) shape": {194: 2, "198a": 1, 217: 1},
    "free 9..15, every seed":                  {194: 2, "198a": 1, 217: 1},
}
DIAM3_198A = {9: 2, 10: 2, 11: 6, 12: 6, 13: 8, 14: 13}

#: 217, exhaustively over degree sequences: how many graphical sequences with
#: max degree <= 6 have residue 2, and the smallest Ls a degree-preserving
#: swap search could reach on any of them. Ls <= 6 is what the hypothesis
#: needs; from n = 11 on, nothing came close.
SEQ_217 = {9: (139, 5), 10: (114, 6), 11: (70, 8), 12: (29, 8),
           13: (7, 10), 14: (1, 11)}


def report(nmax=15, quick=False):
    """How hard the hunt looked, and how close it got."""
    lo, hi = 9, (11 if quick else nmax)
    iters = 400 if quick else 1200
    print("=" * 74)
    print("THE HUNT: minimum defect reached, 9..%d vertices" % hi)
    print("=" * 74)
    print("""
The defect is the hypothesis violated, cleared of its denominator:

    194   n*alpha - n - sum_v l(v)
    198a  n*b - 2n - sum_v ecc(v)
    217   Ls - 4*[residue == 2] - 2

each <= 0 exactly on an applicable graph. A counterexample is a connected
graph with defect <= 0 and no Hamiltonian path.
""")
    print(f"   {'n':>3} {'conj':>6} {'family':>8} {'free':>8}   best structure")
    for n in range(lo, hi + 1):
        sd = seeds(n)
        # block shapes: all singletons for every cut size, and -- the shapes
        # that actually win -- one cut vertex with two singleton blocks and
        # one big one, which is the clique-plus-two-pendants skeleton
        shapes = [(k, [1] * (n - k)) for k in range(1, (n - 2) // 2 + 1)]
        shapes += [(1, [n - 3, 1, 1]), (1, [n - 4, 2, 1]), (1, [n - 5, 2, 2]),
                   (2, [n - 6, 2, 1, 1]), (2, [n - 5, 1, 1, 1])]
        shapes = [(k, s) for k, s in shapes
                  if len(s) >= k + 2 and all(x >= 1 for x in s)]
        for num in (194, "198a", 217):
            fam = min(family_search(num, k, list(s), iters=iters,
                                    seed=n * 31 + k)[0] for k, s in shapes)
            best, tag = 10 ** 9, None
            for name, g0 in sd.items():
                d, _ = free_search(num, g0, iters=iters, seed=n)
                if d < best:
                    best, tag = d, name
            flag = "  *** COUNTEREXAMPLE ***" if min(fam, best) <= 0 else ""
            print(f"   {n:>3} {str(num):>6} {fam:>8} {best:>8}   {tag}{flag}",
                  flush=True)


def verify(g, num):
    """Independent verification of a claimed counterexample: the hypothesis by
    its own definition, and the absence of a Hamiltonian path twice over --
    by the Held-Karp decision and by exhibiting a longest path."""
    n = g.number_of_nodes()
    print(f"   n = {n}, m = {g.number_of_edges()}, connected "
          f"{nx.is_connected(g)}")
    print(f"   graph6 {nx.to_graph6_bytes(g, header=False).strip().decode()}")
    if num == 194:
        print(f"   alpha = {alpha(g)},  1 + l_avg = {1 + avg_l(g):.4f}")
    elif num == "198a":
        print(f"   b = {largest_induced_bipartite(g)},  "
              f"2 + ecc_avg = {2 + avg_ecc(g):.4f}")
    else:
        print(f"   Ls = {max_leaves(g)},  residue = {residue(g)},  "
              f"cap = {4 * (1 if residue(g) == 2 else 0) + 2}")
    print(f"   hypothesis holds: {HYPS[num](g)}")
    p = longest_path(g)
    ok = all(g.has_edge(p[i], p[i + 1]) for i in range(len(p) - 1))
    print(f"   Held-Karp says Hamiltonian path: {has_ham_path(g)}")
    print(f"   longest path has {len(p)} of {n} vertices, genuine path: {ok}")
    print(f"   longest path: {p}")
    return HYPS[num](g) and not has_ham_path(g)


def records():
    """The defect floor, and the graphs that attain it, re-derived here.

    Each of these is non-Hamiltonian and misses its hypothesis by the smallest
    integer margin anything reached at any order. They are what a counter-
    example would have to beat."""
    print("=" * 74)
    print("THE FLOOR: the non-Hamiltonian graphs closest to each hypothesis")
    print("=" * 74)
    rows = []
    for a in (3, 4, 5, 6):
        rows.append((f"K_{{{a},{a + 2}}}", nx.complete_bipartite_graph(a, a + 2),
                     194))
    for m in (6, 8, 10, 12):
        g = nx.complete_graph(m)
        g.add_edges_from([(0, m), (0, m + 1)])
        rows.append((f"K_{m} + 2 pendants at one vertex", g, "198a"))
    for legs in (3, 5, 7):
        g = nx.Graph()
        for i in range(3):
            prev = 0
            for j in range(legs):
                g.add_edge(prev, 1 + i * legs + j)
                prev = 1 + i * legs + j
        rows.append((f"spider, 3 legs of length {legs}", g, 217))
    print(f"   {'graph':<34} {'n':>3} {'conj':>6} {'defect':>7}  Ham path")
    ok = True
    for name, g, num in rows:
        d = DEFECTS[num](g)
        h = has_ham_path(g)
        ok &= (not h)
        print(f"   {name:<34} {g.number_of_nodes():>3} {str(num):>6} "
              f"{d:>7}  {h}")
    print("\n   the defect of each family is constant in n: 2 for 194, 1 for "
          "198a,\n   1 for 217. The fractional slack of the same graphs is "
          "2/n, 1/n, 1,\n   which is why a slack-based hunt looks like it is "
          "converging and is not.")
    return ok


def main():
    full = len(sys.argv) > 1 and sys.argv[1].isdigit()
    print("WOWII 194, 198a, 217 -- counterexample hunt above the exhaustive "
          "range\n")
    print("validating the bitmask invariants against wowii.py")
    assert validate(), "the fast invariants disagree with wowii.py"

    print("\nthe control: Prajapati's 11-vertex counterexample to conjecture "
          "200\nsatisfies none of these three hypotheses, and not narrowly")
    g = nx.from_graph6_bytes(b"J??FFBRq}N_")
    print(f"   194   alpha {alpha(g)} vs 1 + l_avg {1 + avg_l(g):.3f}"
          f"   defect {defect_194(g)}")
    print(f"   198a  b {largest_induced_bipartite(g)} vs 2 + ecc_avg "
          f"{2 + avg_ecc(g):.3f}   defect {defect_198a(g)}")
    print(f"   217   Ls {max_leaves(g)} vs cap "
          f"{4 * (1 if residue(g) == 2 else 0) + 2}   defect {defect_217(g)}")
    print("   so it is no help as a seed, and the hunt starts from scratch")

    print()
    characterise()
    print()
    assert records(), "a record holder turned out to have a Hamiltonian path"

    print("\n" + "=" * 74)
    print("EXHAUSTIVE RANGE")
    print("=" * 74)
    nmax = int(sys.argv[1]) if full else 8
    tot, app, nonham, mind, ce = exhaustive(nmax)
    print(f"\n   {tot} connected graphs on 2..{nmax} vertices")
    for num in (194, "198a", 217):
        d, wit = mind[num]
        print(f"   conj {str(num):>5}  applicable {app[num]:>7}   "
              f"minimum defect over the {nonham} non-Hamiltonian graphs: {d}"
              f"   ({nx.to_graph6_bytes(wit, header=False).strip().decode()})")
    assert not ce, f"COUNTEREXAMPLES FOUND: {ce}"
    print(f"   counterexamples: {len(ce)}")
    if not full:
        r = EXHAUSTIVE[9]
        print(f"\n   and, recorded from the 9-vertex run "
              f"(`python wowii_ham_hunt.py 9` reproduces it):")
        print(f"   {r['graphs']} connected graphs, {r['nonham']} of them "
              f"non-Hamiltonian, {r['counterexamples']} counterexamples")
        print(f"   applicable {r['applicable']}, minimum defect "
              f"{r['min_defect']}")

    print()
    report(quick=not full)

    print("\n" + "=" * 74)
    print("""VERDICT

No counterexample. That is a NULL RESULT and it is weak: conjecture 200
survived a stronger null result than this one -- exhaustive to 8 vertices and
fully covered by five classical certificates -- and is false at 11.

What is not weak:

  * 217 is TRUE for every graph on 15 or more vertices, because its
    hypothesis is unsatisfiable there. See characterise(). Only 10 <= n <= 14
    is open, and over every one of the degree sequences that range permits
    (%s in total) a degree-preserving search never got Ls below 6 for
    n >= 11 at all.

  * 198a has no counterexample of diameter <= 2, proved in characterise().
    Forcing diameter >= 3 makes the search WORSE as n grows -- defect
    %s -- because b(G) grows faster than the average
    eccentricity does once the graph has to branch.

  * The three defect floors are 2, 1 and 1, and they are FLAT in n: the same
    three families (K_{a,a+2}, a clique with two pendants at one vertex, a
    three-legged spider) attain them at every order from 4 to 16. Nothing in
    any search ever got below them.""" % (
        sum(v[0] for v in SEQ_217.values()),
        ", ".join(f"{k}:{v}" for k, v in sorted(DIAM3_198A.items()))))


if __name__ == "__main__":
    main()
