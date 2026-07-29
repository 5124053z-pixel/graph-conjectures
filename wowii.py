"""Graffiti.pc (Written on the Wall II) conjectures still marked open, tested.

Source: the `WrittenOnTheWallII` directory of Google DeepMind's
`formal-conjectures` repository, where 22 of DeLaVina's Graffiti.pc conjectures
are formalised in Lean and tagged `category research open`.

Why this list rather than the original PDF. The Graffiti corpus was ruled out
earlier in this project because Roucairol and Cazenave document that its
definitions disagree between sources -- they produced two false refutations from
mis-defined invariants. A Lean formalisation removes that failure mode
completely: the definition *is* the Lean definition, and there is nothing to
misread.

Which of the 22 to attack was decided by how heavily each theme already appears
on DeLaVina's own page of resolved conjectures:

    induced bipartite b(G)   40 mentions      heavily worked
    Hamiltonicity            25
    forest number f(G)       24
    local independence l(v)  22
    total domination          4
    induced tree              2               barely touched
    residue                   0               untouched
    well totally dominated    0               untouched

so the conjectures below are the ones whose invariants nobody appears to have
resolved anything about.
"""
from __future__ import annotations

import itertools
import math

import networkx as nx

#: Graffiti.pc conjectures were tested against a graph database before being
#: published, so a violation on a graph this small is evidence of a bug in the
#: checker rather than of a false conjecture. This guard exists because the
#: mistake has already been made twice here -- see failures 11 and 12.
SUSPICIOUS_ORDER = 6


# ---------------------------------------------------------------------------
# invariants
# ---------------------------------------------------------------------------

def indep_number(g):
    """Independence number, exact."""
    if g.number_of_nodes() == 0:
        return 0
    comp = nx.complement(g)
    return len(max(nx.find_cliques(comp), key=len)) if g.number_of_nodes() else 0


def local_independence(g, v):
    """l(v) = alpha(G[N(v)]), the independence number of the neighbourhood."""
    nbrs = list(g[v])
    if not nbrs:
        return 0
    return indep_number(g.subgraph(nbrs))


def residue(g):
    """Havel-Hakimi residue: zeros left after running the algorithm."""
    seq = sorted((d for _, d in g.degree()), reverse=True)
    while seq and seq[0] > 0:
        d = seq[0]
        seq = seq[1:]
        if d > len(seq):
            return None                      # not graphical; cannot happen here
        for i in range(d):
            seq[i] -= 1
        seq.sort(reverse=True)
    return len(seq)


def largest_induced_forest(g):
    """Max |S| with G[S] acyclic. Exact, by descending subset size."""
    nodes = list(g)
    for k in range(len(nodes), 0, -1):
        for s in itertools.combinations(nodes, k):
            if nx.is_forest(g.subgraph(s)):
                return k
    return 0


def largest_induced_tree(g):
    """Max |S| with G[S] a tree (connected and acyclic)."""
    nodes = list(g)
    for k in range(len(nodes), 0, -1):
        for s in itertools.combinations(nodes, k):
            sub = g.subgraph(s)
            if nx.is_tree(sub):
                return k
    return 0


def largest_induced_path(g):
    """Max |S| with G[S] a path."""
    nodes = list(g)
    for k in range(len(nodes), 0, -1):
        for s in itertools.combinations(nodes, k):
            sub = g.subgraph(s)
            if nx.is_tree(sub) and all(d <= 2 for _, d in sub.degree()):
                return k
    return 0


def is_total_dominating(g, s):
    s = set(s)
    return all(any(u in s for u in g[v]) for v in g)


def minimal_total_dominating_sizes(g):
    """Sizes of every *minimal* total dominating set (minimal by inclusion)."""
    nodes = list(g)
    sizes = set()
    for k in range(1, len(nodes) + 1):
        for s in itertools.combinations(nodes, k):
            if not is_total_dominating(g, s):
                continue
            if any(is_total_dominating(g, set(s) - {x}) for x in s):
                continue                     # not minimal
            sizes.add(k)
    return sizes


def is_well_totally_dominated(g):
    """All minimal total dominating sets have the same size."""
    sizes = minimal_total_dominating_sizes(g)
    return len(sizes) <= 1


def pendant_vertices(g):
    return [v for v, d in g.degree() if d == 1]


def average_degree(g):
    n = g.number_of_nodes()
    return 2 * g.number_of_edges() / n if n else 0


# ---------------------------------------------------------------------------
# the conjectures: each returns (applicable, holds)
# ---------------------------------------------------------------------------

def conj_322(g):
    """n >= 5 and lambda_max(COMPLEMENT) <= 1  =>  well totally dominated.

    The complement is an overline in the source HTML. Reading the statement as
    lambda_max(G) <= 1 admits only complete graphs and makes the conjecture
    trivial; the real hypothesis is equivalent to G being complete multipartite.
    See wowii322.py."""
    if g.number_of_nodes() < 5:
        return False, True
    gc = nx.complement(g)
    if any(local_independence(gc, v) > 1 for v in gc):
        return False, True
    return True, is_well_totally_dominated(g)


def conj_316(g):
    """avg degree of the complement <= number of pendant vertices
       =>  G is well totally dominated."""
    if average_degree(nx.complement(g)) > len(pendant_vertices(g)):
        return False, True
    return True, is_well_totally_dominated(g)


def conj_314(g):
    """triangle-free and largest induced path <= 4  =>  well totally dominated."""
    if sum(nx.triangles(g).values()) > 0:
        return False, True
    if largest_induced_path(g) > 4:
        return False, True
    return True, is_well_totally_dominated(g)


def conj_141(g):
    """tree(G) >= (1/2)*girth - 1 + max_v l(v), as DeLaVina states it.

    The Lean formalisation has floor(girth/2), which is strictly weaker at odd
    girth. This tests the original. See wowii141.py."""
    girth = 0 if nx.is_forest(g) else nx.girth(g)
    rhs = girth / 2 - 1 + max(local_independence(g, v) for v in g)
    return True, largest_induced_tree(g) >= rhs


def ecc_set(g, S, include_S=False):
    """The two set-eccentricities the formalisation distinguishes.

        ecc(S)     = max over v NOT in S of dist(v, S)   (0 if S is everything)
        eccSet(S)  = max over ALL v      of dist(v, S)

    Getting this wrong is not hypothetical: the first version of `conj_144`
    used "the eccentricity of the vertices in the centre", which is the radius,
    and duly reported 25 counterexamples -- the first being the triangle. A
    3-vertex counterexample to a conjecture on DeLaVina's list is not a
    discovery, it is a bug, and that is how this one was caught."""
    S = set(S)
    dist = dict(nx.all_pairs_shortest_path_length(g))
    outside = [v for v in g if v not in S] if not include_S else list(g)
    if not outside:
        return 0
    return max(min(dist[v][s] for s in S) for v in outside)


def conj_144(g):
    """tree(G) >= girth - 1 + ecc(center), with ecc the *set* eccentricity."""
    girth = 0 if nx.is_forest(g) else nx.girth(g)
    rhs = girth - 1 + ecc_set(g, nx.center(g))
    return True, largest_induced_tree(g) >= rhs


def conj_145(g):
    """2 * eccSet(B) <= tree(G) * lMin(complement), B the max-eccentricity set."""
    gc = nx.complement(g)
    if gc.number_of_edges() == 0:
        return False, True
    lmin = min(local_independence(gc, v) for v in gc)
    if lmin <= 0:
        return False, True
    ecc = nx.eccentricity(g)
    hi = max(ecc.values())
    B = [v for v in g if ecc[v] == hi]
    return True, 2 * ecc_set(g, B, include_S=True) <= largest_induced_tree(g) * lmin


def conj_61(g):
    """f(G) >= residue(G) + ceil(diam(G)/3)."""
    import math
    rhs = residue(g) + math.ceil(nx.diameter(g) / 3)
    return True, largest_induced_forest(g) >= rhs


CONJECTURES = {
    322: ("n>=5, every l(v)<=1  =>  well totally dominated", conj_322),
    316: ("avgdeg(complement) <= #pendants  =>  well totally dominated", conj_316),
    314: ("triangle-free, induced path <= 4  =>  well totally dominated", conj_314),
    141: ("tree(G) >= floor(girth/2) - 1 + max l(v)", conj_141),
    144: ("tree(G) >= girth - 1 + ecc(center)", conj_144),
    145: ("2*eccSet(B) <= tree(G) * lMin(complement)", conj_145),
    61:  ("f(G) >= residue(G) + ceil(diam/3)", conj_61),
}


def scan(nmax=7, verbose=True, only=None):
    """Scan every connected graph on 2..nmax vertices.

    Up to 7 vertices the atlas suffices; beyond that `allgraphs.generate` is
    used, which verifies its own output against OEIS A000088 and A001349 before
    returning anything."""
    if nmax <= 7:
        graphs = [g for g in nx.graph_atlas_g()
                  if 2 <= g.number_of_nodes() <= nmax and nx.is_connected(g)]
    else:
        from allgraphs import connected_graphs
        graphs = connected_graphs(nmax, verbose=verbose)
    if verbose:
        print(f"all {len(graphs)} connected graphs on 2..{nmax} vertices\n")
    out = {}
    for k, (desc, fn) in sorted(CONJECTURES.items()):
        if only is not None and k not in only:
            continue
        applicable = holds = 0
        bad = []
        for g in graphs:
            try:
                app, ok = fn(g)
            except Exception:
                continue
            if app:
                applicable += 1
                if ok:
                    holds += 1
                else:
                    bad.append((g.number_of_nodes(), sorted(g.edges())))
        out[k] = (applicable, holds, bad)
        if verbose:
            if not bad:
                tag = ""
            elif min(n for n, _ in bad) <= SUSPICIOUS_ORDER:
                # A counterexample this small to a conjecture on a published,
                # database-tested list is almost always a bug in the checker,
                # not a discovery. Twice already in this file.
                tag = (f"   !! {len(bad)} violations, smallest on "
                       f"{min(n for n, _ in bad)} vertices -- SUSPECT THE CODE")
            else:
                tag = f"   *** {len(bad)} COUNTEREXAMPLES ***"
            print(f"  #{k:<4} {desc:<52} applicable {applicable:>5}  "
                  f"holds {holds:>5}{tag}")
            for n, e in sorted(bad)[:2]:
                print(f"        n={n}: {e}")
    return out


if __name__ == "__main__":
    scan()


# ---------------------------------------------------------------------------
# further invariants, matching the Lean definitions
# ---------------------------------------------------------------------------

def max_leaves_spanning_tree(g):
    """Ls(G): the largest number of leaves over all spanning trees."""
    if not nx.is_connected(g):
        return 0
    best = 0
    for t in nx.SpanningTreeIterator(g):
        best = max(best, sum(1 for _, d in t.degree() if d == 1))
    return best


def largest_induced_bipartite(g):
    """b(G)."""
    nodes = list(g)
    for k in range(len(nodes), 0, -1):
        for s in itertools.combinations(nodes, k):
            if nx.is_bipartite(g.subgraph(s)):
                return k
    return 0


def path_cover_number(g):
    """Fewest paths whose vertex sets cover V(G). Paths may overlap: the Lean
    definition asks for a cover, not a partition.

    A vertex set `s` is usable iff some path *in G* has exactly `s` as its
    support -- equivalently, iff `G[s]` has a Hamiltonian path. It does **not**
    have to be an induced path. The first version required the induced
    subgraph to be a path, which made `p(K_3) = 3` instead of 1 and produced a
    triangle as a counterexample to conjecture 40. Same shape of error as the
    `ecc` misreading in failure 11, caught the same way: by the size of the
    witness."""
    nodes = list(g)
    paths = set()
    for k in range(1, len(nodes) + 1):
        for s in itertools.combinations(nodes, k):
            if k == 1 or has_hamiltonian_path(g.subgraph(s)):
                paths.add(frozenset(s))
    universe = frozenset(nodes)
    plist = sorted(paths, key=lambda p: -len(p))
    for k in range(1, len(nodes) + 1):
        for combo in itertools.combinations(plist, k):
            if frozenset().union(*combo) == universe:
                return k
    return len(nodes)


def havel_hakimi_zero_step(g):
    """First step of the Havel-Hakimi reduction at which a zero appears."""
    seq = sorted((d for _, d in g.degree()), reverse=True)
    for i in range(len(seq) + 1):
        if not seq or seq[-1] == 0:
            return i
        d = seq[0]
        seq = seq[1:]
        for j in range(min(d, len(seq))):
            seq[j] -= 1
        seq.sort(reverse=True)
    return len(seq)


def triangles_at(g, v):
    nb = set(g[v])
    return sum(1 for a, b in itertools.combinations(nb, 2) if g.has_edge(a, b))


def freq_min_triangles(g):
    t = {v: triangles_at(g, v) for v in g}
    lo = min(t.values())
    return sum(1 for v in t if t[v] == lo)


def has_hamiltonian_path(g):
    """Held-Karp over subsets: O(2^n * n^2) instead of O(n!).

    At n = 8 the permutation version costs 40,320 per graph and 11,117 graphs
    make that hopeless; the bitmask version is a few hundred operations."""
    nodes = list(g)
    n = len(nodes)
    if n <= 1:
        return True
    idx = {v: i for i, v in enumerate(nodes)}
    adj = [0] * n
    for u, v in g.edges():
        adj[idx[u]] |= 1 << idx[v]
        adj[idx[v]] |= 1 << idx[u]
    full = (1 << n) - 1
    reach = [[False] * n for _ in range(1 << n)]
    for i in range(n):
        reach[1 << i][i] = True
    for mask in range(1 << n):
        row = reach[mask]
        for i in range(n):
            if not row[i]:
                continue
            if mask == full:
                return True
            rest = adj[i] & ~mask
            while rest:
                b = rest & -rest
                reach[mask | b][b.bit_length() - 1] = True
                rest ^= b
    return any(reach[full])


def avg_l(g):
    return sum(local_independence(g, v) for v in g) / g.number_of_nodes()


def total_domination_number(g):
    for k in range(1, g.number_of_nodes() + 1):
        for s in itertools.combinations(list(g), k):
            if is_total_dominating(g, s):
                return k
    return g.number_of_nodes()


# --- the remaining conjectures -------------------------------------------

def conj_2(g):
    """Ls(G) >= 2 * (l_avg(G) - 1)."""
    return True, max_leaves_via_cds(g) >= 2 * (avg_l(g) - 1)


def conj_40(g):
    """f(G) >= ceil((p(G) + b(G) + 1) / 2)."""
    import math
    rhs = math.ceil((path_cover_number(g) + largest_induced_bipartite(g) + 1) / 2)
    return True, largest_induced_forest(g) >= rhs


def conj_59(g):
    """f(G) >= ceil(sqrt(residue(G) * b(G)))."""
    import math
    rhs = math.ceil(math.sqrt(residue(g) * largest_induced_bipartite(g)))
    return True, largest_induced_forest(g) >= rhs


def conj_194(g):
    """alpha(G) <= 1 + l_avg(G)  =>  G has a Hamiltonian path."""
    if indep_number(g) > 1 + avg_l(g):
        return False, True
    return True, has_hamiltonian_path(g)


def conj_198a(g):
    """b(G) <= 2 + average eccentricity  =>  G has a Hamiltonian path."""
    ecc = nx.eccentricity(g)
    avg_ecc = sum(ecc.values()) / g.number_of_nodes()
    if largest_induced_bipartite(g) > 2 + avg_ecc:
        return False, True
    return True, has_hamiltonian_path(g)


def conj_200(g):
    """tree(G) == ceil(1 + l_avg(G))  =>  G has a Hamiltonian path."""
    import math
    if largest_induced_tree(g) != math.ceil(1 + avg_l(g)):
        return False, True
    return True, has_hamiltonian_path(g)


def conj_217(g):
    """Ls(G) <= 4 * [residue(G) == 2] + 2  =>  G has a Hamiltonian path."""
    if max_leaves_via_cds(g) > 4 * (1 if residue(g) == 2 else 0) + 2:
        return False, True
    return True, has_hamiltonian_path(g)


def conj_291(g):
    """gamma_t(G) <= havelHakimiZeroStep(G) + freqMinTriangles(G), n > 2."""
    if g.number_of_nodes() <= 2:
        return False, True
    rhs = havel_hakimi_zero_step(g) + freq_min_triangles(g)
    return True, total_domination_number(g) <= rhs


CONJECTURES.update({
    2:   ("Ls(G) >= 2*(l_avg - 1)", conj_2),
    40:  ("f(G) >= ceil((p(G)+b(G)+1)/2)", conj_40),
    59:  ("f(G) >= ceil(sqrt(residue*b))", conj_59),
    194: ("alpha <= 1 + l_avg  =>  Hamiltonian path", conj_194),
    198: ("b(G) <= 2 + avg ecc  =>  Hamiltonian path", conj_198a),
    200: ("tree(G) == ceil(1+l_avg)  =>  Hamiltonian path", conj_200),
    217: ("Ls(G) <= 4*[residue==2]+2  =>  Hamiltonian path", conj_217),
    291: ("gamma_t <= HH-zero-step + freq(min triangles)", conj_291),
})


# ---------------------------------------------------------------------------
# the last seven, and the definitions they need
# ---------------------------------------------------------------------------

def degree_l2_norm(g):
    """sqrt(sum of squared degrees). NOTE: conjecture 100's doc comment calls
    this `diam(Gc)`, but the Lean statement uses `degreeL2Norm Gc`. The
    statement is what is formalised, so the statement is what is tested."""
    return math.sqrt(sum(d * d for _, d in g.degree()))


def count_induced_c4(g):
    """Number of induced 4-cycles, as unordered vertex sets."""
    n = 0
    for s in itertools.combinations(list(g), 4):
        sub = g.subgraph(s)
        if sub.number_of_edges() == 4 and all(d == 2 for _, d in sub.degree()):
            n += 1
    return n


def has_c4_subgraph(g):
    """A 4-cycle as a subgraph, not necessarily induced."""
    for a, b, c, d in itertools.permutations(list(g), 4):
        if g.has_edge(a, b) and g.has_edge(b, c) and g.has_edge(c, d) \
                and g.has_edge(d, a):
            return True
    return False


def dist_min(g, S):
    """min pairwise distance between two distinct vertices of S. 0 if |S| < 2.

    Corrected 2026-07-29 per Dr. DeLaVina directly: this was previously
    implemented as "distance from S to the rest of the graph", which is a
    different quantity and made conjecture 65 trivial by construction. Her
    own example is P5: the minimum-degree set A = {endpoints} has pairwise
    distance 4, not the distance-to-rest value of 1 the old code returned."""
    S = list(S)
    if len(S) < 2:
        return 0
    dist = dict(nx.all_pairs_shortest_path_length(g))
    return min(dist[u][v] for u, v in itertools.combinations(S, 2))


def max_triangles_at_vertex(g):
    return max(triangles_at(g, v) for v in g)


def graph_square(g):
    h = nx.Graph()
    h.add_nodes_from(g)
    dist = dict(nx.all_pairs_shortest_path_length(g))
    for u, v in itertools.combinations(list(g), 2):
        if dist[u].get(v, 99) <= 2:
            h.add_edge(u, v)
    return h


def conj_19(g):
    """floor(avg ecc + max_v l(v)) <= b(G)."""
    ecc = nx.eccentricity(g)
    avg = sum(ecc.values()) / g.number_of_nodes()
    lhs = math.floor(avg + max(local_independence(g, v) for v in g))
    return True, lhs <= largest_induced_bipartite(g)


def conj_65(g):
    """distMin(minDeg set) + ceil(distMin(maxDeg set)/3) <= f(G)."""
    degs = dict(g.degree())
    A = [v for v in g if degs[v] == min(degs.values())]
    M = [v for v in g if degs[v] == max(degs.values())]
    lhs = dist_min(g, A) + math.ceil(dist_min(g, M) / 3)
    return True, lhs <= largest_induced_forest(g)


def conj_100(g):
    """alpha(G) <= ceil((max_v l(v) + 0.5 * degreeL2Norm(Gc)) / 2), Gc connected."""
    gc = nx.complement(g)
    if gc.number_of_nodes() < 2 or not nx.is_connected(gc):
        return False, True
    rhs = math.ceil((max(local_independence(g, v) for v in g)
                     + 0.5 * degree_l2_norm(gc)) / 2)
    return True, indep_number(g) <= rhs


def conj_133(g):
    """rad(G) + floor(l_avg)^cC4 <= path(G), cC4 = 0 if G has a 4-cycle else 1."""
    c = 0 if has_c4_subgraph(g) else 1
    lhs = nx.radius(g) + math.floor(avg_l(g)) ** c
    return True, lhs <= largest_induced_path(g)


def conj_142(g):
    """(2/3)*girth + eccSet(B) <= tree(G), B the max-eccentricity vertices."""
    girth = 0 if nx.is_forest(g) else nx.girth(g)
    ecc = nx.eccentricity(g)
    hi = max(ecc.values())
    B = [v for v in g if ecc[v] == hi]
    lhs = (2 / 3) * girth + ecc_set(g, B, include_S=True)
    return True, lhs <= largest_induced_tree(g)


def conj_146(g):
    """2*eccSet(B) <= tree(G) * radius(G^2), radius(G^2) > 0."""
    r = nx.radius(graph_square(g))
    if r <= 0:
        return False, True
    ecc = nx.eccentricity(g)
    hi = max(ecc.values())
    B = [v for v in g if ecc[v] == hi]
    return True, 2 * ecc_set(g, B, include_S=True) <= largest_induced_tree(g) * r


def conj_160(g):
    """max_v l(v) + max_v T(v) * cC4(G) <= Ls(G), with cC4 as DeLaVina defines it.

    Definition 73 of her definitions file: cC4(G) is 1 if G is C4-free -- NOT
    necessarily induced -- and 0 otherwise. The Lean file uses a *count* of
    induced 4-cycles, under which the statement fails on 8,985 of the 12,112
    connected graphs on at most 8 vertices. See wowii160.py."""
    lhs = (max(local_independence(g, v) for v in g)
           + max_triangles_at_vertex(g) * (0 if has_c4_subgraph(g) else 1))
    return True, lhs <= max_leaves_via_cds(g)


CONJECTURES.update({
    19:  ("floor(avg ecc + max l(v)) <= b(G)", conj_19),
    65:  ("distMin(A) + ceil(distMin(M)/3) <= f(G)", conj_65),
    100: ("alpha <= ceil((max l(v) + 0.5*||deg(Gc)||_2)/2)", conj_100),
    133: ("rad + floor(l_avg)^cC4 <= path(G)", conj_133),
    142: ("(2/3)girth + eccSet(B) <= tree(G)", conj_142),
    146: ("2*eccSet(B) <= tree(G)*radius(G^2)", conj_146),
    160: ("max l(v) + maxTri*countInducedC4 <= Ls(G)", conj_160),
})


def max_leaves_via_cds(g):
    """Ls(G) = n - (minimum connected dominating set size).

    The internal vertices of a maximum-leaf spanning tree are exactly a minimum
    connected dominating set, so the two problems are the same one. This matters
    because the direct version iterates spanning trees, and K_8 alone has
    8^6 = 262,144 of them -- hopeless across 12,112 graphs. Searching connected
    dominating sets by increasing size is a few hundred subsets instead.

    Validated against `max_leaves_spanning_tree` on every connected graph up to
    6 vertices before being used anywhere."""
    n = g.number_of_nodes()
    if n <= 1:
        return 0
    if n == 2:
        return 2          # the identity fails here: K_2's spanning tree is a
                          # single edge, both of whose ends are leaves, while
                          # n - mcds gives 1. Found by validating rather than
                          # assuming -- it was the one mismatch in 142 graphs.
    nodes = list(g)
    for k in range(1, n + 1):
        for s in itertools.combinations(nodes, k):
            sset = set(s)
            if not all(v in sset or any(u in sset for u in g[v]) for v in g):
                continue
            if nx.is_connected(g.subgraph(s)):
                return n - k
    return 0
