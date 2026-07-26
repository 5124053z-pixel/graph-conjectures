"""Counterexample search over equitable partitions, for bounds 44 and 46.

Bounds 44 and 46 are the only two candidates in the family whose status is
unknown (Damnjanovic-Ha-Stevanovic, arXiv:2606.14550, June 2026). Everything
else is settled.

The method is the one that refuted the other twelve, and it rests on a
simplification worth stating plainly:

    **the entire counterexample condition depends only on the quotient matrix.**

If a graph G has an equitable partition into cells C_1..C_k with quotient matrix
B = (b_ij) -- every vertex of C_i has exactly b_ij neighbours in C_j -- then every
vertex of C_i has the same degree s_i = sum_j b_ij and the same neighbour-average
m_i = (sum_j b_ij s_j) / s_i. So the right-hand side of the conjecture is a
function of B alone. And mu(G) >= rho(L_B) for L_B = diag(s) - B (Lemma 2.2
there, standard). Hence

    rho(L_B) > max_{b_ij > 0} f(s_i, s_j, m_i, m_j)

certifies a counterexample, with no graph ever constructed and no bound on its
size. Cell sizes enter only through realizability.

Realizability (their Lemma 2.3): cell sizes n_1..n_k work iff b_ii <= n_i - 1
with b_ii or n_i even, b_ij <= n_j, and n_i b_ij = n_j b_ji. Multiplying every
n_i by a common even factor preserves the last condition and makes the first
three easier, so **realizability reduces to the ratio system n_i/n_j = b_ji/b_ij
being consistent** -- and then the graph exists, however large it has to be.

Parametrising by cell sizes makes consistency automatic: fix n, and for i < j put
b_ij = c * n_j / gcd(n_i, n_j), b_ji = c * n_i / gcd(n_i, n_j) for an integer
c in [0, gcd(n_i, n_j)]. Every consistent B arises this way.
"""
from __future__ import annotations

import argparse
import math
import random
from math import gcd

import numpy as np

NEG_INF = float("-inf")

#: Slack below this counts as a counterexample; anything in [TOL, 0] is exact
#: equality seen through floating point.
#:
#: This family is *exactly tight* on every regular bipartite graph, so a naive
#: `slack < 0` test reports K_{3,3}, C_4, Q_3 and infinitely many others as
#: counterexamples at the 1e-16 level. That already happened twice in this
#: repository -- see failures 2 and 7 in the README -- which is why the number
#: lives here once instead of being retyped in each scan.
TOL = -1e-9


def is_counterexample(slack_value):
    return slack_value is not None and slack_value < TOL


# ---------------------------------------------------------------------------
# the two open bounds, with the paper's convention: a negative radicand makes
# the term -inf, i.e. the edge drops out of the maximum. Clamping to 0 instead
# would inflate the bound and hide counterexamples.
# ---------------------------------------------------------------------------

def bound_44(a, b, p, q):
    r = 2 * ((a - 1) ** 2 + (b - 1) ** 2 + p * q - a * b)
    return 2 + math.sqrt(r) if r >= 0 else NEG_INF


def bound_46(a, b, p, q):
    r = 2 * (a * a + b * b) - 16 * a * b / (p + q) + 4
    return 2 + math.sqrt(r) if r >= 0 else NEG_INF


BOUNDS = {44: bound_44, 46: bound_46}


# ---------------------------------------------------------------------------
# quotient matrices
# ---------------------------------------------------------------------------

def degrees(B):
    """(s_i, m_i) for each cell: degree, and average degree of neighbours."""
    s = B.sum(axis=1)
    if (s <= 0).any():
        return None, None
    m = (B @ s) / s
    return s, m


def quotient_rhs(B, f):
    """The conjecture's right-hand side for any graph with this quotient."""
    s, m = degrees(B)
    if s is None:
        return None
    k = len(s)
    best = NEG_INF
    for i in range(k):
        for j in range(i, k):
            if B[i, j] > 0:
                v = f(s[i], s[j], m[i], m[j])
                if v > best:
                    best = v
    return best


def mu_lower(B):
    """rho(L_B) <= mu(G) for any graph with this quotient matrix."""
    s = B.sum(axis=1)
    return float(np.abs(np.linalg.eigvals(np.diag(s) - B)).max())


def connected(B):
    """Is the quotient graph connected? Cells are adjacent when b_ij > 0."""
    k = len(B)
    seen, stack = {0}, [0]
    while stack:
        i = stack.pop()
        for j in range(k):
            if j not in seen and (B[i, j] > 0 or B[j, i] > 0):
                seen.add(j)
                stack.append(j)
    return len(seen) == k


def slack(B, f, prune=True):
    """rhs - mu_lower. Negative certifies a counterexample.

    Two prunings, both of which remove only non-counterexamples but between them
    remove the plateau that otherwise traps the search at slack exactly 0:

    * **regular quotients.** If every cell has the same degree s the graph is
      s-regular, every bound in the family equals 2s, and mu <= 2s. Slack is
      identically 0 there.
    * **disconnected quotients.** The right-hand side is a maximum over all
      edges and mu is a maximum over components, so a single regular bipartite
      component pins the slack at exactly 0 no matter what the rest looks like.
      This was what the first version of the search kept rediscovering. Any
      counterexample can be taken connected (pass to the component of largest
      mu), so requiring connectivity loses nothing."""
    s, _ = degrees(B)
    if s is None:
        return None
    if prune and (s.max() - s.min() < 1e-12 or not connected(B)):
        return None
    rhs = quotient_rhs(B, f)
    if rhs is None:
        return None
    return rhs - mu_lower(B)


def build(n, c):
    """Quotient matrix from cell sizes n and coefficients c[(i,j)], i<j, plus
    diagonal entries c[(i,i)]. Guarantees n_i b_ij = n_j b_ji by construction."""
    k = len(n)
    B = np.zeros((k, k))
    for i in range(k):
        B[i, i] = c.get((i, i), 0)
        for j in range(i + 1, k):
            g = gcd(n[i], n[j])
            cij = c.get((i, j), 0)
            B[i, j] = cij * n[j] // g
            B[j, i] = cij * n[i] // g
    return B


def realizable(n, B):
    """The remaining conditions of Lemma 2.3, after scaling n by an even factor
    (which is always allowed). Returns the scale factor, or None."""
    k = len(n)
    t = 2
    for _ in range(64):
        ok = True
        for i in range(k):
            if B[i, i] > t * n[i] - 1:
                ok = False
            for j in range(k):
                if i != j and B[i, j] > t * n[j]:
                    ok = False
        if ok:
            return t
        t *= 2
    return None


# ---------------------------------------------------------------------------
# validation against the twelve published counterexamples
# ---------------------------------------------------------------------------

#: bound -> (cell sizes, quotient matrix), Table 3 of arXiv:2606.14550
PUBLISHED = {
    11: ((1, 7, 21), [[0, 7, 0], [1, 0, 3], [0, 1, 2]]),
    13: ((1, 9, 36), [[0, 9, 0], [1, 0, 4], [0, 1, 3]]),
    18: ((48, 12, 42, 56), [[0, 1, 0, 42], [4, 0, 14, 0],
                            [0, 4, 0, 4], [36, 0, 3, 0]]),
    19: ((5, 5, 5, 5), [[0, 5, 1, 0], [5, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]),
    20: ((5, 5, 5, 5), [[0, 5, 1, 0], [5, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]]),
    21: ((69, 19, 76), [[0, 0, 76], [0, 0, 8], [69, 2, 0]]),
    22: ((7, 9, 14), [[0, 9, 2], [7, 0, 0], [1, 0, 0]]),
    24: ((69, 19, 76), [[0, 0, 76], [0, 0, 8], [69, 2, 0]]),
    30: ((83, 38, 44, 88), [[0, 0, 0, 88], [0, 0, 22, 0],
                            [0, 19, 0, 6], [83, 0, 3, 0]]),
    40: ((1, 39, 312), [[0, 39, 0], [1, 0, 8], [0, 1, 0]]),
    47: ((28, 2, 22, 7), [[0, 2, 0, 1], [28, 0, 0, 0],
                          [0, 0, 0, 7], [4, 0, 22, 0]]),
    56: ((1, 5, 15, 15), [[0, 5, 0, 0], [1, 0, 3, 0], [0, 1, 0, 1], [0, 0, 1, 0]]),
}

from laplacian_bounds import TYPE1, TYPE2  # noqa: E402


def _published_f(k):
    """Wrap a bound from the transcribed list with the -inf convention."""
    if k in TYPE1:
        g = TYPE1[k]
        return lambda a, b, p, q: max(_neg(g, a, p), _neg(g, b, q))
    g = TYPE2[k]
    return lambda a, b, p, q: _neg(g, a, b, p, q)


def _neg(g, *args):
    try:
        v = g(*args)
    except (ValueError, ZeroDivisionError):
        return NEG_INF
    return v if isinstance(v, float) and math.isfinite(v) else NEG_INF


def validate(verbose=True):
    """Every published counterexample must come out as a counterexample here.
    This checks the quotient machinery, the -inf convention and the bound
    transcriptions all at once, against twelve independent known answers."""
    bad = []
    for k, (n, rows) in sorted(PUBLISHED.items()):
        B = np.array(rows, dtype=float)
        f = _published_f(k)
        s, m = degrees(B)
        rhs, mu = quotient_rhs(B, f), mu_lower(B)
        t = realizable(n, B)
        ok = rhs is not None and rhs - mu < TOL and t is not None
        if not ok:
            bad.append((k, rhs, mu, t))
        if verbose:
            print(f"bound {k:>2}: rhs = {rhs:>10.6f}   mu >= {mu:>10.6f}   "
                  f"slack = {rhs - mu:>+10.6f}   {'OK' if ok else 'FAILED'}")
    return bad


# ---------------------------------------------------------------------------
# search
# ---------------------------------------------------------------------------

def random_state(k, nmax, rng):
    n = [rng.randint(1, nmax) for _ in range(k)]
    c = {}
    for i in range(k):
        c[(i, i)] = rng.choice([0, 0, 0, rng.randint(0, 4)])
        for j in range(i + 1, k):
            g = gcd(n[i], n[j])
            c[(i, j)] = rng.randint(0, g) if rng.random() < 0.6 else 0
    return n, c


def perturb(n, c, k, nmax, rng):
    n, c = list(n), dict(c)
    what = rng.random()
    if what < 0.45:
        i = rng.randrange(k)
        n[i] = max(1, min(nmax, n[i] + rng.choice([-3, -2, -1, 1, 2, 3])))
    elif what < 0.9:
        i = rng.randrange(k)
        j = rng.randrange(k)
        i, j = min(i, j), max(i, j)
        lim = (n[i] - 1) if i == j else gcd(n[i], n[j])
        c[(i, j)] = max(0, min(lim, c.get((i, j), 0)
                               + rng.choice([-2, -1, 1, 2])))
    else:
        i, j = rng.randrange(k), rng.randrange(k)
        i, j = min(i, j), max(i, j)
        c[(i, j)] = 0 if c.get((i, j), 0) else 1
    for i in range(k):
        for j in range(i, k):
            lim = (n[i] - 1) if i == j else gcd(n[i], n[j])
            if c.get((i, j), 0) > lim:
                c[(i, j)] = lim
    return n, c


def anneal(conj, k, nmax=400, iters=60_000, restarts=40, seed=0, verbose=True):
    """Minimise slack = rhs - mu over realizable quotient matrices."""
    f = BOUNDS[conj]
    rng = random.Random(seed)
    best = (math.inf, None, None)

    for r in range(restarts):
        n, c = random_state(k, nmax, rng)
        cur = slack(build(n, c), f)
        if cur is None:
            cur = math.inf
        temp0 = 2.0
        for it in range(iters):
            T = temp0 * (1 - it / iters) + 1e-3
            n2, c2 = perturb(n, c, k, nmax, rng)
            B2 = build(n2, c2)
            s2 = slack(B2, f)
            if s2 is None:
                continue
            if s2 < cur or rng.random() < math.exp(-(s2 - cur) / T):
                n, c, cur = n2, c2, s2
                if cur < best[0] and realizable(n2, B2) is not None:
                    best = (cur, list(n2), B2.copy())
        if verbose and (r + 1) % 10 == 0:
            print(f"    restart {r + 1:>3}/{restarts}   best slack = {best[0]:+.6f}")
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--conjecture", type=int, default=None, choices=[44, 46])
    ap.add_argument("--cells", type=int, default=4)
    ap.add_argument("--nmax", type=int, default=400)
    ap.add_argument("--iters", type=int, default=60_000)
    ap.add_argument("--restarts", type=int, default=40)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--validate-only", action="store_true")
    args = ap.parse_args()

    print("validating the quotient machinery against the twelve published "
          "counterexamples\n")
    bad = validate()
    print()
    if bad:
        print(f"VALIDATION FAILED for bounds {[b[0] for b in bad]} -- not "
              f"searching until this is fixed")
        raise SystemExit(1)
    print("all twelve reproduced\n")
    if args.validate_only:
        return

    targets = [args.conjecture] if args.conjecture else [44, 46]
    for conj in targets:
        print(f"searching for a counterexample to bound {conj} "
              f"({args.cells} cells, n <= {args.nmax})")
        s, n, B = anneal(conj, args.cells, args.nmax, args.iters,
                         args.restarts, args.seed)
        print(f"  best slack (rhs - mu) = {s:+.12f}")
        if TOL <= s <= -TOL:
            print("  (exactly tight, not a counterexample)")
        if is_counterexample(s):
            print(f"  COUNTEREXAMPLE   cell sizes {n}")
            print(f"  quotient matrix\n{B.astype(int)}")
        else:
            print("  no counterexample found (not a proof)")
        print()


if __name__ == "__main__":
    main()
