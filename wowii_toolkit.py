"""The free lower bounds, and how far they take each Graffiti.pc conjecture.

Every partial result in this project came from one move: find a bound on the
invariant that costs nothing to prove, then see which range of the conjecture it
already covers. This file collects the bounds, each with its one-line reason,
and reports the residual each conjecture is left with.

None of these is deep. What they are is *unused*: the conjectures had been
sitting on DeLaVina's list since Graffiti.pc produced them, and nobody had
checked which of them a star, a shortest path or a chordless cycle already
settles.
"""
from __future__ import annotations

import itertools
import math

import networkx as nx

from domination import domination_number as DOM

import wowii as W


def girth(g):
    return 0 if nx.is_forest(g) else nx.girth(g)


def maxl(g):
    return max(W.local_independence(g, v) for v in g)


# --- the bounds -----------------------------------------------------------

def b_acyclic(g):
    """A forest is its own largest induced tree/forest."""
    return g.number_of_nodes() if nx.is_forest(g) else 0


def b_star(g):
    """v plus a maximum independent set of N(v) is a star."""
    return 1 + maxl(g)


def b_path(g):
    """A shortest path is induced."""
    return nx.diameter(g) + 1


def b_cycle(g):
    """A shortest cycle is chordless; drop one vertex."""
    return max(girth(g) - 1, 0)


def shortest_cycles(g):
    gi = girth(g)
    if gi == 0:
        return []
    out = []
    for s in itertools.combinations(list(g), gi):
        sub = g.subgraph(s)
        if sub.number_of_edges() == gi and all(d == 2 for _, d in sub.degree()) \
                and nx.is_connected(sub):
            out.append(set(s))
    return out


def b_caterpillar_cycle(g):
    """C - x is an induced path; hang an independent set attaching once each."""
    gi = girth(g)
    if gi < 4:
        return 0
    best = 0
    for C in shortest_cycles(g):
        for x in C:
            P = C - {x}
            cand = [a for a in g if a not in C and len(set(g[a]) & P) == 1]
            m = W.indep_number(g.subgraph(cand)) if cand else 0
            best = max(best, (gi - 1) + m)
    return best


def b_caterpillar_indep(g):
    """I independent, plus an independent set attaching once each: a star forest."""
    comp = nx.complement(g)
    best = 0
    for I in nx.find_cliques(comp):
        Is = set(I)
        cand = [a for a in g if a not in Is and len(set(g[a]) & Is) == 1]
        m = W.indep_number(g.subgraph(cand)) if cand else 0
        best = max(best, len(I) + m)
    return best


def b_two_indep(g):
    """Two disjoint independent sets induce a bipartite subgraph."""
    comp = nx.complement(g)
    best = 0
    for I in nx.find_cliques(comp):
        rest = [v for v in g if v not in set(I)]
        j = W.indep_number(g.subgraph(rest)) if rest else 0
        best = max(best, len(I) + j)
    return best


def b_cyclepath(g):
    """girth >= 4: hang a shortest path off a shortest cycle.

    Complementary to b_caterpillar_cycle, not subsumed by it: that one attaches
    single vertices at several cycle vertices, this one attaches a whole path at
    one. Dropping it when the toolkit was consolidated is what left conjecture
    144 with a residual."""
    gi = girth(g)
    if gi < 4:
        return 0
    cs = shortest_cycles(g)
    if not cs:
        return 0
    dist = dict(nx.all_pairs_shortest_path_length(g))
    best = 0
    for C in cs:
        e = max(min(dist[v][x] for x in C) for v in g)
        best = max(best, gi - 1 + e)
    return best


def b_one_vertex_off(g):
    """If deleting a single vertex leaves a tree, then tree(G) >= n - 1.

    Not circular: the test is one connectivity-and-acyclicity check per vertex,
    O(n*m), while computing tree(G) is exponential. It just names a witness that
    happens to be easy to look for."""
    n = g.number_of_nodes()
    for v in g:
        h = nx.restricted_view(g, [v], [])
        if nx.is_connected(h) and nx.is_forest(h):
            return n - 1
    return 0


def b_rad_tree(g):
    """tree(G) >= 2*rad - 1.

    Erdos, Saks and Sos, *Maximum induced trees in graphs*, JCTB 41 (1986),
    Theorem 2.2 -- stated there for the induced path number, with a proof due
    to Fan Chung, and an induced path is an induced tree."""
    return 2 * nx.radius(g) - 1


def b_three_step_star(g):
    """girth >= 6: the star at v extends two steps without closing a cycle."""
    return 3 + maxl(g) if girth(g) >= 6 else 0


def b_two_step_star(g):
    """girth >= 5: the star at v extends one step further without closing a cycle."""
    return 2 + maxl(g) if girth(g) >= 5 else 0


def b_cyclomatic(g):
    """Delete one vertex per independent cycle: m - n + 1 of them suffice."""
    return g.number_of_nodes() - (g.number_of_edges() - g.number_of_nodes() + 1)


def b_rad_path(g):
    """path(G) >= 2*rad - 1.

    Erdos, Saks and Sos, *Maximum induced trees in graphs*, JCTB 41 (1986),
    Theorem 2.2, with a proof due to Fan Chung. Not a free observation --
    a published theorem, found only when the bounds were checked against
    the literature rather than only against the corpus."""
    return 2 * nx.radius(g) - 1


def b_half_bipartite(g):
    """f(G) >= b(G)/2 + 1 for connected G on more than one vertex.

    Not found here: DeLaVina states it in the note under conjecture 40 on her
    own page, dated 6 March 2004, together with the corollary that 40 follows
    whenever the path covering number is 1. Verified with zero violations."""
    return W.largest_induced_bipartite(g) / 2 + 1


# --- witness certificates --------------------------------------------------
# These differ in kind from everything above. A structural bound (the star, the
# ball, Gallai) is an ingredient of a proof for *all* graphs. A witness only
# exhibits one object and reports its size, which certifies the conjecture on
# the graph in hand and says nothing about any other graph. They finish the
# in-range verification; they are not steps toward a theorem.
def b_greedy_cds(g):
    """Ls >= n - |S| for any connected dominating set S we can exhibit."""
    n = g.number_of_nodes()
    best = n
    for start in g:
        S = {start}
        dom = set(g[start]) | {start}
        while dom != set(g):
            cand = [u for v in S for u in g[v] if u not in S]
            u = max(cand, key=lambda u: len((set(g[u]) | {u}) - dom))
            S.add(u)
            dom |= set(g[u]) | {u}
        best = min(best, len(S))
    return n - best


def b_greedy_path(g):
    """path(G) >= |P| for any induced path P we can exhibit."""
    best = 1
    for s in g:
        P, used = [s], {s}
        for end in (-1, 0):
            while True:
                ext = [u for u in g[P[end]] if u not in used
                       and sum(1 for x in P if g.has_edge(u, x)) == 1]
                if not ext:
                    break
                u = min(ext, key=g.degree)
                P.append(u) if end == -1 else P.insert(0, u)
                used.add(u)
        best = max(best, len(P))
    return best


def b_exhibit_forest(g):
    """f(G) >= |I u A| for I a maximum independent set and A independent with
    the bipartite graph between them acyclic. Generalises the star-forest
    bound, whose condition -- one neighbour in I each -- is sufficient for
    acyclicity but not necessary."""
    best = 0
    nodes = list(g)
    for I in nx.find_cliques(nx.complement(g)):
        rest = [v for v in nodes if v not in I]
        for k in range(len(rest), 0, -1):
            if len(I) + k <= best:
                break
            for A in itertools.combinations(rest, k):
                if any(g.has_edge(x, y) for x, y in itertools.combinations(A, 2)):
                    continue
                if nx.is_forest(g.subgraph(list(I) + list(A))):
                    best = max(best, len(I) + k)
                    break
    return best


TREE_BOUNDS = [b_acyclic, b_star, b_path, b_cycle, b_caterpillar_cycle,
               b_cyclepath, b_two_step_star, b_three_step_star,
               b_one_vertex_off, b_rad_tree]
FOREST_BOUNDS = [b_acyclic, b_path, b_caterpillar_indep, b_cyclomatic,
                 lambda g: W.indep_number(g) + 1,
                 W.largest_induced_tree]
FOREST_BOUNDS = FOREST_BOUNDS + [b_one_vertex_off, b_rad_tree,
                                 b_exhibit_forest, b_half_bipartite]
BIPARTITE_BOUNDS = FOREST_BOUNDS + [b_two_indep, W.largest_induced_forest]
PATH_BOUNDS = [b_path, b_cycle, b_rad_path, b_greedy_path]


# --- the conjectures ------------------------------------------------------

def T_eccB_marker(g):
    return eccB(g)


def eccB(g):
    e = nx.eccentricity(g)
    hi = max(e.values())
    return W.ecc_set(g, [v for v in g if e[v] == hi], include_S=True)


def avg_ecc(g):
    e = nx.eccentricity(g)
    return sum(e.values()) / g.number_of_nodes()


# --- bounds for the upper-bound conjectures --------------------------------
def b_maxdeg_leaves(g):
    """Ls(G) >= Delta: root a spanning tree at a maximum-degree vertex."""
    return max(d for _, d in g.degree())


def u_gallai(g):
    """alpha <= n - nu (Gallai): a maximum matching's edges each keep a vertex
    out of any independent set."""
    return g.number_of_nodes() - len(nx.max_weight_matching(g))


def u_tdom_maxdeg(g):
    """gamma_t <= n - Delta + 1 for connected G on >= 3 vertices."""
    return g.number_of_nodes() - max(d for _, d in g.degree()) + 1


def u_tdom_edge(g):
    """gamma_t <= 2 when some edge uv has N(u) u N(v) = V -- an exact witness,
    found in O(m*n) instead of solving the set cover."""
    for u, v in g.edges():
        if set(g[u]) | set(g[v]) == set(g):
            return 2
    return g.number_of_nodes()


def lmin_complement(g):
    gc = nx.complement(g)
    return min(W.local_independence(gc, v) for v in gc)


def b_tree_over_lmin(g):
    """For 145: `2*eccB <= tree * lMin(Gc)`.

    If lMin(Gc) >= 2 the right side is at least 2*tree >= 2*(diam+1), and
    eccB <= diam always, so it holds outright. Otherwise lMin = 1 and the
    statement is just `2*eccB <= tree`, which the tree toolkit answers."""
    lm = lmin_complement(g)
    if lm >= 2:
        return 2 * (nx.diameter(g) + 1) * lm
    return max(b(g) for b in TREE_BOUNDS)


def u_tdom_twice(g):
    """gamma_t <= 2*gamma: each dominator drags in one neighbour."""
    return 2 * DOM(g)


CONJ = {
    "19":  (BIPARTITE_BOUNDS, lambda g: math.floor(avg_ecc(g) + maxl(g))),
    "40":  (FOREST_BOUNDS,    lambda g: math.ceil((W.path_cover_number(g)
                                                   + W.largest_induced_bipartite(g) + 1) / 2)),
    "59":  (FOREST_BOUNDS,    lambda g: math.ceil(math.sqrt(
                                          W.residue(g) * W.largest_induced_bipartite(g)))),
    "61":  (FOREST_BOUNDS,    lambda g: W.residue(g) + math.ceil(nx.diameter(g) / 3)),
    "133": (PATH_BOUNDS,      lambda g: nx.radius(g) + math.floor(W.avg_l(g))
                                        ** (0 if W.has_c4_subgraph(g) else 1)),
    "141": (TREE_BOUNDS,      lambda g: girth(g) // 2 - 1 + maxl(g)),
    "142": (TREE_BOUNDS,      lambda g: (2 / 3) * girth(g) + eccB(g)),
    "144": (TREE_BOUNDS,      lambda g: girth(g) - 1 + W.ecc_set(g, nx.center(g))),
    "146": (TREE_BOUNDS,      lambda g: 2 * eccB(g)
                                        / max(nx.radius(W.graph_square(g)), 1)),
    "2":   ([b_maxdeg_leaves, b_greedy_cds], lambda g: 2 * (W.avg_l(g) - 1)),
    "145": ([b_tree_over_lmin], lambda g: 2 * T_eccB_marker(g)),
}

# Conjectures stated as `lhs <= rhs`: prove them by squeezing an upper bound on
# the left through. Same method, mirrored.
CONJ_UPPER = {
    "100": ([u_gallai],
            lambda g: math.ceil((maxl(g)
                                 + 0.5 * W.degree_l2_norm(nx.complement(g))) / 2),
            lambda g: W.conj_100(g)[0]),
    "291": ([u_tdom_maxdeg, u_tdom_twice, u_tdom_edge],
            lambda g: W.havel_hakimi_zero_step(g) + W.freq_min_triangles(g),
            lambda g: W.conj_291(g)[0]),
}


def main():
    graphs = [g for g in nx.graph_atlas_g()
              if 2 <= g.number_of_nodes() <= 7 and nx.is_connected(g)]
    print(f"free-bound coverage over the {len(graphs)} connected graphs on "
          f"2..7 vertices\n")
    print(f"   {'conj':>5} {'covered':>9} {'residual':>9}")
    for name, (bounds, rhs) in sorted(CONJ.items(), key=lambda t: int(t[0])):
        left = [g for g in graphs if not any(b(g) >= rhs(g) for b in bounds)]
        print(f"   {name:>5} {len(graphs) - len(left):>9} {len(left):>9}"
              f"{'   <-- settled in range' if not left else ''}")
    for name, (bounds, rhs, applies) in sorted(CONJ_UPPER.items(),
                                               key=lambda t: int(t[0])):
        ap = [g for g in graphs if applies(g)]
        left = [g for g in ap if not any(b(g) <= rhs(g) for b in bounds)]
        print(f"   {name:>5} {len(ap) - len(left):>9} {len(left):>9}"
              f"   (of {len(ap)} where the hypothesis applies)"
              f"{'   <-- settled in range' if not left else ''}")


if __name__ == "__main__":
    main()
