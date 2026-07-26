"""Classical, *proved* upper bounds on the Laplacian spectral radius.

The point of this file is the other half of the problem. The literature attacks
these 68 conjectures exclusively by search -- deep cross-entropy, Monte Carlo
tree search, exhaustive generation -- i.e. everyone is trying to *refute* them.
Nobody appears to have checked the cheap direction: whether a conjectured bound
is simply implied by a theorem that has been known since the 1980s.

If a conjectured bound B satisfies B(G) >= M(G) for every graph G, where M is a
proved upper bound on mu, then the conjecture is true and needs no search at
all. Since all these expressions are local, the implication can usually be
established pointwise, in one line of algebra.

Each bound below is stated with its source and then *numerically re-validated*
against a corpus of graphs by `validate()`, on the same principle as the rest of
this repository: a misremembered theorem would otherwise produce a confident
false proof, which is worse than no result.
"""
from __future__ import annotations

import math

import networkx as nx

from laplacian_bounds import laplacian_spectral_radius, local_degrees

SQ = math.sqrt


# --- type 1: max over vertices of f(d_v, m_v) ------------------------------

KNOWN1 = {
    # Merris, "A note on Laplacian graph eigenvalues", LAA 285 (1998) 33-35.
    "merris": lambda d, m: d + m,
    # Li & Pan, "De Caen's inequality and bounds on the largest Laplacian
    # eigenvalue of a graph", LAA 328 (2001) 153-160.
    "li-pan": lambda d, m: SQ(2 * d * d + 2 * d * m),
    # Zhang, "Two sharp upper bounds for the Laplacian eigenvalues",
    # LAA 376 (2004) 207-213.
    "zhang": lambda d, m: d + SQ(d * m),
}

# --- type 2: max over edges of f(d_u, d_v, m_u, m_v) -----------------------

KNOWN2 = {
    # Anderson & Morley, "Eigenvalues of the Laplacian of a graph",
    # Linear and Multilinear Algebra 18 (1985) 141-145.
    "anderson-morley": lambda a, b, p, q: a + b,
    # Das, "An improved upper bound for Laplacian graph eigenvalues", LAA 368
    # (2003) 269-278.
    "das": lambda a, b, p, q: (a * (a + p) + b * (b + q)) / (a + b),
    # The Perron-vector bound extracted from the proof of Proposition 3.5 in
    # Damnjanovic-Ha-Stevanovic (arXiv:2606.14550). They use it only to confirm
    # bound 42; it is stated there as a step, not as a bound, but it is one.
    #
    # With Q = D + A the signless Laplacian and P = diag(d_i), the matrix
    # B = P^{-1} Q P is similar to Q. Take a positive Perron vector v of B, let
    # i maximise v_i and let j be the neighbour of i with the largest v_j. The
    # Perron equations give (rho - d_i) v_i <= m_i v_j and (rho - d_j) v_j <=
    # m_j v_i, and rho is at least every diagonal entry, so
    # (rho - d_i)(rho - d_j) <= m_i m_j. Solving the quadratic and using
    # mu <= rho gives the bound below at that one edge, hence at the maximum
    # over all edges.
    "perron": lambda a, b, p, q: (a + b + SQ((a - b) ** 2 + 4 * p * q)) / 2,
}


def evaluate(name, g):
    d, m = local_degrees(g)
    if name in KNOWN1:
        f = KNOWN1[name]
        return max(f(d[v], m[v]) for v in g)
    f = KNOWN2[name]
    return max(f(d[u], d[v], m[u], m[v]) for u, v in g.edges())


def corpus(max_atlas=7, max_tree=14, randoms=400, seed=1):
    """A corpus wide enough that a misremembered theorem fails on it."""
    import random
    for g in nx.graph_atlas_g():
        if len(g) >= 2 and nx.is_connected(g):
            yield g
    for n in range(3, max_tree + 1):
        for g in nx.nonisomorphic_trees(n):
            yield g
    rng = random.Random(seed)
    for _ in range(randoms):
        n = rng.randint(5, 22)
        g = nx.gnp_random_graph(n, rng.uniform(0.08, 0.9), seed=rng.randint(0, 1 << 30))
        if nx.is_connected(g):
            yield g
    # deliberately irregular shapes: stars, double stars, brooms, kites
    for n in range(4, 40):
        yield nx.star_graph(n - 1)
        yield nx.lollipop_graph(max(3, n // 2), n - max(3, n // 2))
        for a in range(1, n - 2):
            g = nx.Graph()
            g.add_edges_from([(0, 1)])
            g.add_edges_from((0, 2 + i) for i in range(a))
            g.add_edges_from((1, 2 + a + i) for i in range(n - 2 - a))
            if nx.is_connected(g):
                yield g


def validate(tol=1e-9, verbose=True):
    """Check every claimed theorem against the corpus. Returns {name: worst
    violation}; an entry means the bound is NOT valid as transcribed and must
    not be used for proofs."""
    names = list(KNOWN1) + list(KNOWN2)
    worst = {n: 0.0 for n in names}
    witness = {}
    count = 0
    for g in corpus():
        if g.number_of_edges() == 0:
            continue
        count += 1
        mu = laplacian_spectral_radius(g)
        for n in names:
            v = evaluate(n, g)
            if mu - v > worst[n]:
                worst[n] = mu - v
                witness[n] = sorted(g.edges())
    bad = {n: (w, witness[n]) for n, w in worst.items() if w > tol}
    if verbose:
        print(f"validated {len(names)} classical bounds on {count:,} graphs")
        for n in names:
            status = "OK" if worst[n] <= tol else f"INVALID (mu exceeds by {worst[n]:.6g})"
            print(f"  {n:<18} {status}")
    return bad


if __name__ == "__main__":
    bad = validate()
    raise SystemExit(1 if bad else 0)
