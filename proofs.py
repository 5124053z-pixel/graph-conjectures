"""Conjectures resolved as TRUE, with the pointwise inequality each proof rests on.

Every proof here has the same one-line shape:

    conjectured bound  >=  a classical proved bound   (pointwise in d and m)

so the conjecture follows immediately. `verify()` re-checks each pointwise
inequality numerically over the feasible region, so that a typo in a proof
shows up as a test failure rather than as a wrong claim in the README.

The proofs themselves are in the docstrings; they are short enough to check by
eye, which is the point -- these are conjectures that survived deep
cross-entropy search and eight Monte Carlo search algorithms, not because they
are hard, but because everyone was searching for counterexamples to statements
that are provable in one line.
"""
from __future__ import annotations

import math

SQ = math.sqrt

#: Conjectures settled by the max-degree argument. Let v* have maximum degree
#: Delta. Every neighbour of v* has degree <= Delta, so m_{v*} <= Delta =
#: d_{v*}: at a vertex of maximum degree, d >= m always. Anderson-Morley gives
#: mu <= max_{uv in E}(d_u + d_v) <= 2*Delta. So it suffices that the bound is
#: >= 2d on the half-plane m <= d -- much weaker than pointwise domination,
#: since it only has to hold where m <= d.
#:
#: conjecture -> (bound, the inequality it reduces to on m <= d)
MAXDEG = {
    1:  (lambda d, m: SQ(4 * d**3 / m), "4d^3/m >= 4d^2  <=>  d >= m"),
    4:  (lambda d, m: 2 * d * d / m,    "2d^2/m >= 2d  <=>  d >= m"),
    7:  (lambda d, m: d * d / m + d,    "d^2/m >= d  <=>  d >= m"),
    14: (lambda d, m: 2 * d**3 / m**2,  "d^3/m^2 >= d  <=>  d^2 >= m^2  <=>  d >= m"),
    16: (lambda d, m: 2 * d**4 / m**3,  "d^4/m^3 >= d  <=>  d^3 >= m^3  <=>  d >= m"),
    5:  (lambda d, m: d * d / m + m,    "d^2/m + m >= 2d by AM-GM, for every d, m"),
}

#: conjecture -> (theorem used, pointwise inequality as a function that must be
#: >= 0 everywhere on the feasible region, human-readable justification)
PROVED = {
    5: (
        "anderson-morley",
        lambda d, m: d * d / m + m - 2 * d,
        """d^2/m + m - 2d = (d - m)^2 / m >= 0, so the bound is >= 2 d_v at every
        vertex, hence >= 2*Delta at a vertex of maximum degree. Anderson-Morley
        gives mu <= max_{uv in E} (d_u + d_v) <= 2*Delta.""",
    ),
    6: (
        "li-pan",
        lambda d, m: (m * m + 3 * d * d) - 2 * d * (d + m),
        """(m^2 + 3d^2) - 2d(d + m) = (d - m)^2 >= 0, so
        sqrt(m^2 + 3d^2) >= sqrt(2 d (d + m)), and Li-Pan gives
        mu <= max_v sqrt(2 d_v^2 + 2 d_v m_v).""",
    ),
    9: (
        "li-pan",
        lambda d, m: ((m + 3 * d) / 2) ** 2 - 2 * d * (d + m),
        """((m + 3d)/2)^2 - 2d(d + m) = (d - m)^2 / 4 >= 0, so
        (m + 3d)/2 >= sqrt(2 d (d + m)); Li-Pan.""",
    ),
    12: (
        "merris",
        lambda d, m: (2 * m * m + 2 * d * d) - (d + m) ** 2,
        """2m^2 + 2d^2 - (d + m)^2 = (d - m)^2 >= 0, so
        sqrt(2m^2 + 2d^2) >= d + m; Merris gives mu <= max_v (d_v + m_v).""",
    ),
    34: (
        "anderson-morley",
        lambda a, b: 2 * (a * a + b * b) - (a + b) ** 2,
        """2(a^2 + b^2) - (a + b)^2 = (a - b)^2 >= 0, so
        2(d_u^2 + d_v^2)/(d_u + d_v) >= d_u + d_v; Anderson-Morley.""",
    ),
    37: (
        "anderson-morley",
        lambda a, b: 2 * (a * a + b * b) - (a + b) ** 2,
        """Same inequality: sqrt(2(d_u^2 + d_v^2)) >= d_u + d_v;
        Anderson-Morley.""",
    ),
    38: (
        "anderson-morley",
        lambda a, b: 2 * ((a - 1) ** 2 + (b - 1) ** 2) - ((a - 1) + (b - 1)) ** 2,
        """With x = d_u - 1 >= 0 and y = d_v - 1 >= 0, sqrt(2(x^2 + y^2)) >= x + y,
        so 2 + sqrt(2(d_u-1)^2 + 2(d_v-1)^2) >= 2 + (d_u - 1) + (d_v - 1)
        = d_u + d_v; Anderson-Morley.""",
    ),
}


def verify(steps=4001, hi=80.0, tol=-1e-9, verbose=True):
    """Re-check every pointwise inequality. Returns the list of failures."""
    bad = []
    for k, (thm, ineq, _) in sorted(PROVED.items()):
        worst, at = math.inf, None
        for i in range(1, steps):
            x = hi * i / steps
            for j in range(1, steps, 7):
                y = hi * j / steps
                v = ineq(x, y)
                if v < worst:
                    worst, at = v, (x, y)
        if worst < tol:
            bad.append((k, thm, worst, at))
        if verbose:
            print(f"conjecture {k:>2}  via {thm:<16} "
                  f"min of pointwise inequality = {worst:+.3e}"
                  f"{'   FAIL' if worst < tol else ''}")
    return bad


def verify_maxdeg(steps=200_001, tol=-1e-9, verbose=True):
    """Each MAXDEG bound must satisfy B(d, m) >= 2d on m <= d. The expressions
    are homogeneous of degree 1, so setting d = 1 and sweeping m in (0, 1] is
    exhaustive."""
    bad = []
    for k, (f, why) in sorted(MAXDEG.items()):
        worst, at = math.inf, None
        for i in range(1, steps + 1):
            m = i / steps
            v = f(1.0, m) - 2.0
            if v < worst:
                worst, at = v, m
        if worst < tol:
            bad.append((k, worst, at))
        if verbose:
            print(f"conjecture {k:>2}  via mu <= 2*Delta   "
                  f"min(B - 2d) on m <= d = {worst:+.3e}"
                  f"{'   FAIL' if worst < tol else ''}")
    return bad


def resolved():
    return sorted(set(PROVED) | set(MAXDEG))


if __name__ == "__main__":
    n = len(resolved())
    print(f"{n} conjectures proved true\n")
    print("-- by domination of a classical bound --")
    bad = verify()
    print()
    print("-- by the max-degree argument (mu <= 2*Delta) --")
    bad += verify_maxdeg()
    print()
    print("all inequalities hold" if not bad else f"FAILURES: {bad}")
    print(f"resolved: {resolved()}")
    raise SystemExit(1 if bad else 0)
