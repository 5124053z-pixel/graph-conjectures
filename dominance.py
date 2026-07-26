"""Prove open conjectures by showing they are implied by a classical theorem.

If B is a conjectured upper bound on mu and M is a *proved* one, then
B(G) >= M(G) for all G implies the conjecture. Because every expression here is
local, that reduces to a pointwise inequality between two small algebraic
functions, which can be scanned numerically and then verified by hand.

The scan is deliberately conservative: it searches for a *violation* of the
pointwise inequality over the feasible parameter region, and only reports a
conjecture as a proof candidate when no violation is found. A candidate still
has to be checked algebraically before it is claimed -- see PROVED in
`results.py` for the ones that survived that step.

Feasible region. For an edge ij write a = d_i, b = d_j, p = m_i, q = m_j. Then
a, b are positive integers and

    a*p = b + (sum of the degrees of the other a-1 neighbours of i)  >= b + a-1

so p >= (a + b - 1)/a, and symmetrically q >= (a + b - 1)/b. When a = 1 the
only neighbour of i is j, so p = b exactly. Using these constraints (rather
than p, q ranging freely) is what makes several of the implications go through.
"""
from __future__ import annotations

import argparse
import math

from laplacian_bounds import TYPE1, TYPE2, OPEN
from known_bounds import KNOWN1, KNOWN2


def _safe(f, *args):
    try:
        v = f(*args)
    except (ZeroDivisionError, ValueError, OverflowError):
        return None
    if isinstance(v, complex) or not math.isfinite(v):
        return None
    return v


# ---------------------------------------------------------------------------
# type 1 conjecture vs type 1 theorem: both homogeneous of degree 1, so the
# inequality depends only on t = d/m. A dense 1-D scan settles it.
# ---------------------------------------------------------------------------

def scan_1v1(fb, fm, steps=200_001, tmax=60.0):
    worst, at = math.inf, None
    for i in range(1, steps):
        t = tmax * i / steps
        b, m = _safe(fb, t, 1.0), _safe(fm, t, 1.0)
        if b is None or m is None:
            continue
        if b - m < worst:
            worst, at = b - m, t
    return worst, at


# ---------------------------------------------------------------------------
# everything else: scan the 4-dimensional edge-parameter region.
# ---------------------------------------------------------------------------

def edge_params(dmax=26, msteps=14):
    """Feasible (a, b, p, q) tuples, respecting the constraints above."""
    for a in range(1, dmax + 1):
        for b in range(1, dmax + 1):
            plo, qlo = (a + b - 1) / a, (a + b - 1) / b
            phi = plo if a == 1 else max(plo, dmax)
            qhi = qlo if b == 1 else max(qlo, dmax)
            ps = [b] if a == 1 else [plo + (phi - plo) * i / msteps
                                     for i in range(msteps + 1)]
            qs = [a] if b == 1 else [qlo + (qhi - qlo) * i / msteps
                                     for i in range(msteps + 1)]
            for p in ps:
                for q in qs:
                    yield a, b, p, q


def scan_4d(value_b, value_m, dmax=26, msteps=14):
    worst, at = math.inf, None
    for a, b, p, q in edge_params(dmax, msteps):
        vb, vm = value_b(a, b, p, q), value_m(a, b, p, q)
        if vb is None or vm is None:
            continue
        if vb - vm < worst:
            worst, at = vb - vm, (a, b, p, q)
    return worst, at


def check(conj, thm, dmax=26, msteps=14):
    """Smallest value of (conjectured bound - proved bound) over the feasible
    region. >= 0 means the conjecture follows from the theorem."""
    in1, kn1 = conj in TYPE1, thm in KNOWN1

    if in1 and kn1:
        return scan_1v1(TYPE1[conj], KNOWN1[thm])

    if in1:                      # type-1 conjecture vs edge-indexed theorem:
        f, m = TYPE1[conj], KNOWN2[thm]      # the max over the two endpoints
        def vb(a, b, p, q):                  # has to beat the theorem's value
            xs = [_safe(f, a, p), _safe(f, b, q)]
            xs = [x for x in xs if x is not None]
            return max(xs) if xs else None
        return scan_4d(vb, lambda a, b, p, q: _safe(m, a, b, p, q), dmax, msteps)

    f = TYPE2[conj]
    if kn1:                      # edge-indexed conjecture vs vertex theorem:
        m = KNOWN1[thm]          # must beat the theorem at the endpoint i
        return scan_4d(lambda a, b, p, q: _safe(f, a, b, p, q),
                       lambda a, b, p, q: _safe(m, a, p), dmax, msteps)

    m = KNOWN2[thm]
    return scan_4d(lambda a, b, p, q: _safe(f, a, b, p, q),
                   lambda a, b, p, q: _safe(m, a, b, p, q), dmax, msteps)


# ---------------------------------------------------------------------------
# A second, sharper criterion: pin the vertex instead of quantifying over all.
#
# The pointwise tests above ask for B_v >= M_v at *every* vertex, which is much
# more than is needed -- only max_v B_v >= mu is required. The gap can be closed
# by naming a vertex where the inequality is easy:
#
#   let v* have maximum degree Delta. Every neighbour of v* has degree <= Delta,
#   so m_{v*} <= Delta = d_{v*}: **at a vertex of maximum degree, d >= m always**.
#   And mu <= 2*Delta (Anderson-Morley). So if
#
#       B(d, m) >= 2d   whenever   m <= d,
#
#   then max_v B_v >= B_{v*} >= 2*Delta >= mu, and the conjecture holds.
#
# This is strictly weaker than pointwise domination -- it only has to hold on
# the half-plane m <= d -- and it settles several conjectures that no theorem
# dominates everywhere.
# ---------------------------------------------------------------------------

def scan_maxdeg_1(f, steps=200_001):
    """min over 0 < m <= d of (f(d, m) - 2d), by scale invariance with d = 1."""
    worst, at = math.inf, None
    for i in range(1, steps + 1):
        m = i / steps                      # m/d in (0, 1]
        v = _safe(f, 1.0, m)
        if v is None:
            continue
        if v - 2.0 < worst:
            worst, at = v - 2.0, m
    return worst, at


def scan_maxdeg_2(f, dmax=26, msteps=16):
    """min over edges incident to a maximum-degree vertex of (f - 2*Delta).

    At such an edge a = Delta, so b <= a, and every m is an average of degrees
    bounded by Delta, so p <= a and q <= a as well."""
    worst, at = math.inf, None
    for a in range(1, dmax + 1):
        for b in range(1, a + 1):
            plo, qlo = (a + b - 1) / a, (a + b - 1) / b
            ps = [b] if a == 1 else [plo + (a - plo) * i / msteps
                                     for i in range(msteps + 1) if plo <= a]
            qs = [a] if b == 1 else [qlo + (a - qlo) * i / msteps
                                     for i in range(msteps + 1) if qlo <= a]
            for p in ps:
                for q in qs:
                    v = _safe(f, a, b, p, q)
                    if v is None:
                        continue
                    if v - 2 * a < worst:
                        worst, at = v - 2 * a, (a, b, p, q)
    return worst, at


def check_maxdeg(conj, dmax=26, msteps=16):
    if conj in TYPE1:
        return scan_maxdeg_1(TYPE1[conj])
    return scan_maxdeg_2(TYPE2[conj], dmax, msteps)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dmax", type=int, default=26)
    ap.add_argument("--msteps", type=int, default=14)
    ap.add_argument("--tol", type=float, default=-1e-9)
    args = ap.parse_args()

    theorems = list(KNOWN1) + list(KNOWN2)
    print(f"testing {len(OPEN)} open conjectures against {len(theorems)} "
          f"classical bounds\n")

    hits = {}
    for k in OPEN:
        row = []
        for t in theorems:
            worst, at = check(k, t, args.dmax, args.msteps)
            row.append((t, worst, at))
            if worst >= args.tol:
                hits.setdefault(k, []).append((t, worst, at))
        best = max(row, key=lambda r: r[1])
        mark = "  <== IMPLIED BY " + best[0].upper() if best[1] >= args.tol else ""
        print(f"conj {k:>2}: closest theorem {best[0]:<16} "
              f"min(B - M) = {best[1]:+.6f}{mark}")

    print()
    print("second criterion: B >= 2d on the half-plane m <= d "
          "(evaluated at a maximum-degree vertex, where mu <= 2*Delta)\n")
    for k in OPEN:
        worst, at = check_maxdeg(k, args.dmax, args.msteps)
        mark = "  <== TRUE" if worst >= args.tol else ""
        if worst >= args.tol:
            hits.setdefault(k, []).append(("mu <= 2*Delta", worst, at))
        print(f"conj {k:>2}: min(B - 2d) over m <= d = {worst:+.6f}{mark}")

    print()
    if hits:
        print(f"{len(hits)} conjecture(s) appear to follow from a known theorem:")
        for k, lst in sorted(hits.items()):
            for t, worst, at in lst:
                print(f"  conjecture {k:>2}  <=  {t}   (margin {worst:+.3g})")
        print("\nThese are CANDIDATES. Each must be verified algebraically "
              "before being claimed.")
    else:
        print("no conjecture is implied pointwise by any transcribed theorem.")


if __name__ == "__main__":
    main()
