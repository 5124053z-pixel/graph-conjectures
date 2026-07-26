"""The 68 conjectured upper bounds on the Laplacian spectral radius.

Source chain:

  * the conjectures were generated automatically by Brankov, Hansen and
    Stevanovic from a database of 273,214 connected graphs on at most 9
    vertices;
  * Al-Yakoob, Ghebleh, Kanso and Stevanovic attacked all 68 with Wagner's deep
    cross-entropy method (arXiv:2509.01607) and refuted 25, plus 5 more by
    exhaustive search over subquartic graphs on <= 14 vertices;
  * Taieb, Roucairol, Cazenave and Harutyunyan (LION 2025) attacked them with
    eight Monte Carlo search algorithms and refuted 29 in total, two of which
    (45 and 48) had been open.

That leaves **36 open**, listed in OPEN below. They all have the form

    mu(G) <= max_v  f(d_v, m_v)          ("type 1")
    mu(G) <= max_{uv in E} f(d_u,d_v,m_u,m_v)   ("type 2")

with d_v the degree of v and m_v the average degree of the neighbours of v.

The formulas here were transcribed from a PDF appendix, so they are checked by
`selftest()` rather than trusted: on a k-regular graph every one of the 68
bounds must evaluate to exactly 2k, and on a k-regular *bipartite* graph
mu = 2k exactly. A transcription error almost surely breaks that identity, so
the self-test is run before any scan.
"""
from __future__ import annotations

import math

import numpy as np
import networkx as nx

SQ = math.sqrt


def local_degrees(g):
    """(d_v, m_v) for every vertex: degree, and average degree of neighbours."""
    d = dict(g.degree())
    m = {v: sum(d[u] for u in g[v]) / d[v] for v in g}
    return d, m


def laplacian_spectral_radius(g):
    a = nx.to_numpy_array(g, nodelist=list(g))
    return float(np.linalg.eigvalsh(np.diag(a.sum(axis=1)) - a)[-1])


# ---------------------------------------------------------------------------
# type 1: max over vertices of f(d, m)
# ---------------------------------------------------------------------------

TYPE1 = {
    1:  lambda d, m: SQ(4 * d**3 / m),
    2:  lambda d, m: 2 * m**2 / d,
    3:  lambda d, m: m**2 / d + m,
    4:  lambda d, m: 2 * d**2 / m,
    5:  lambda d, m: d**2 / m + m,
    6:  lambda d, m: SQ(m**2 + 3 * d**2),
    7:  lambda d, m: d**2 / m + d,
    8:  lambda d, m: SQ(d * (m + 3 * d)),
    9:  lambda d, m: (m + 3 * d) / 2,
    10: lambda d, m: SQ(d * (d + 3 * m)),
    11: lambda d, m: 2 * m**3 / d**2,
    12: lambda d, m: SQ(2 * m**2 + 2 * d**2),
    13: lambda d, m: 2 * m**4 / d**3,
    14: lambda d, m: 2 * d**3 / m**2,
    15: lambda d, m: SQ(4 * m**3 / d),
    16: lambda d, m: 2 * d**4 / m**3,
    17: lambda d, m: (5 * d**4 + 11 * m**4) ** 0.25,
    18: lambda d, m: SQ(2 * m**3 / d + 2 * d**2),
    19: lambda d, m: (4 * d**4 + 12 * d * m**3) ** 0.25,
    20: lambda d, m: SQ(7 * d**2 + 9 * m**2) / 2,
    21: lambda d, m: SQ(d**3 / m + 3 * m**2),
    22: lambda d, m: (2 * d**4 + 14 * d**2 * m**2) ** 0.25,
    23: lambda d, m: SQ(d**2 + 3 * d * m),
    24: lambda d, m: (6 * d**4 + 10 * m**4) ** 0.25,
    25: lambda d, m: (3 * d**4 + 13 * d**2 * m**2) ** 0.25,
    26: lambda d, m: SQ(5 * d**2 + 11 * d * m) / 2,
    27: lambda d, m: SQ((3 * d**2 + 5 * d * m) / 2),
    28: lambda d, m: SQ(2 * m**4 / d**2 + 2 * d * m),
    29: lambda d, m: SQ(m**2 + 3 * m**3 / d),
    30: lambda d, m: m**3 / d**2 + d**2 / m,
    31: lambda d, m: 4 * m**2 / (m + d),
    32: lambda d, m: SQ(m**3 * (m + 3 * d)) / d,
}

# ---------------------------------------------------------------------------
# type 2: max over edges of f(d_i, d_j, m_i, m_j)
# ---------------------------------------------------------------------------


def _sq0(x):
    """sqrt clamped at 0: a few bounds have radicands that go slightly negative
    on very irregular graphs. Clamping makes the bound *smaller*, i.e. easier to
    refute, so it can never manufacture a false positive."""
    return SQ(x) if x > 0 else 0.0


TYPE2 = {
    33: lambda a, b, p, q: 2 * (a + b) - (p + q),
    34: lambda a, b, p, q: 2 * (a**2 + b**2) / (a + b),
    35: lambda a, b, p, q: 2 * (a**2 + b**2) / (p + q),
    36: lambda a, b, p, q: 2 * (p**2 + q**2) / (a + b),
    37: lambda a, b, p, q: SQ(2 * (a**2 + b**2)),
    38: lambda a, b, p, q: 2 + SQ(2 * (a - 1) ** 2 + 2 * (b - 1) ** 2),
    39: lambda a, b, p, q: 2 + _sq0(2 * (a**2 + b**2) - 4 * (p + q) + 4),
    40: lambda a, b, p, q: 2 + _sq0(2 * ((p - 1) ** 2 + (q - 1) ** 2)
                                   + (a**2 + b**2) - (a * p + b * q)),
    41: lambda a, b, p, q: (2 + (p + q) - (a + b)
                            + _sq0(2 * (a**2 + b**2) - 4 * (p + q) + 4)),
    42: lambda a, b, p, q: SQ(a**2 + b**2 + 2 * p * q),
    43: lambda a, b, p, q: 2 + _sq0(3 * (p**2 + q**2) - 2 * p * q
                                    - 4 * (a + b) + 4),
    44: lambda a, b, p, q: 2 + _sq0(2 * ((a - 1) ** 2 + (b - 1) ** 2
                                         + p * q - a * b)),
    45: lambda a, b, p, q: 2 + _sq0((a - b) ** 2 + 2 * (a * p + b * q)
                                    - 4 * (p + q) + 4),
    46: lambda a, b, p, q: 2 + _sq0(2 * (a**2 + b**2)
                                    - 16 * a * b / (p + q) + 4),
    47: lambda a, b, p, q: (2 * (a**2 + b**2) - (p - q) ** 2) / (a + b),
    48: lambda a, b, p, q: (2 * (a**2 + b**2)
                            / (2 + _sq0(2 * (a**2 + b**2) - 4 * (p + q) + 4))),
    49: lambda a, b, p, q: 2 + _sq0(2 * (p**2 + q**2) + (a - b) ** 2
                                    - 4 * (a + b) + 4),
    50: lambda a, b, p, q: 2 * (a**2 + b**2 + p * q - a * b) / (a + b),
    51: lambda a, b, p, q: 2 * (p + q) - 4 * p * q / (a + b),
    52: lambda a, b, p, q: 2 + _sq0(_sq0(8 * (p**4 + q**4)
                                         - 8 * (a**2 + b**2) + 4)
                                    - 4 * (a + b) + 6),
    53: lambda a, b, p, q: 2 + _sq0(_sq0(8 * (p**4 + q**4)
                                         - 8 * (a * p + b * q) + 4)
                                    - 4 * (a + b) + 6),
    54: lambda a, b, p, q: 2 + _sq0(2 * (p**2 + q**2) + (a * p + b * q)
                                    - (a**2 + b**2) - 4 * (a + b) + 4),
    55: lambda a, b, p, q: 2 + _sq0(3 * (p**2 + q**2) - (a**2 + b**2)
                                    - 4 * (p + q) + 4),
    56: lambda a, b, p, q: (a**2 + b**2) * (p + q) / (2 * a * b),
    57: lambda a, b, p, q: 2 + _sq0(2 * (p**2 + q**2)
                                    - 8 * (a**2 + b**2) / (p + q) + 4),
    58: lambda a, b, p, q: 2 + _sq0(2 * (p**2 + p * q + q**2)
                                    - (a * p + b * q) - 4 * (a + b) + 4),
    59: lambda a, b, p, q: (2 * (p**2 + p * q + q**2)
                            - (a**2 + b**2)) / (p + q),
    60: lambda a, b, p, q: 2 + _sq0(2 * (p**2 + p * q + q**2)
                                    - (a**2 + b**2) - 4 * (a + b) + 4),
    61: lambda a, b, p, q: (2 * (p**2 + q**2)
                            / (2 + SQ(2 * ((a - 1) ** 2 + (b - 1) ** 2)))),
    62: lambda a, b, p, q: 2 + _sq0(p**2 + 4 * p * q + q**2 - 2 * a * b
                                    - 4 * (a + b) + 4),
    63: lambda a, b, p, q: a + b + p + q - 4 * a * b / (p + q),
    64: lambda a, b, p, q: p * q * (a + b) / (a * b),
    65: lambda a, b, p, q: (p + q) * (a * p + b * q) / (2 * p * q),
    66: lambda a, b, p, q: (p**2 + 4 * p * q + q**2 - (a * p + b * q)) / (a + b),
    67: lambda a, b, p, q: (p + q) * (a * p + b * q) / (2 * a * b),
    68: lambda a, b, p, q: 2 + _sq0((p - q) ** 2 + 4 * a * b
                                    - 4 * (p + q) + 4),
}

#: still open as of Taieb et al. (LION 2025)
OPEN = [1, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 19, 20, 21, 22, 23,
        24, 25, 26, 27, 30, 33, 34, 35, 37, 38, 39, 40, 42, 44, 46, 47, 56]

#: refuted before this project started
REFUTED = sorted(set(TYPE1) | set(TYPE2) - set()) and [
    k for k in sorted(set(TYPE1) | set(TYPE2)) if k not in OPEN]


def bound(k, g, d=None, m=None):
    """Value of conjectured bound number k on graph g."""
    if d is None:
        d, m = local_degrees(g)
    if k in TYPE1:
        f = TYPE1[k]
        return max(f(d[v], m[v]) for v in g)
    f = TYPE2[k]
    return max(f(d[u], d[v], m[u], m[v]) for u, v in g.edges())


def all_bounds(g, keys=None):
    d, m = local_degrees(g)
    return {k: bound(k, g, d, m) for k in (keys or list(TYPE1) + list(TYPE2))}


def selftest(verbose=False):
    """Every bound must be exactly 2k on a k-regular graph, and mu = 2k on the
    k-regular bipartite ones. Returns a list of discrepancies (empty = pass)."""
    bad = []
    cases = [(nx.complete_bipartite_graph(4, 4), 4, True),
             (nx.cycle_graph(8), 2, True),
             (nx.hypercube_graph(3), 3, True),
             (nx.complete_bipartite_graph(3, 3), 3, True),
             (nx.petersen_graph(), 3, False)]
    for g, k, bipartite in cases:
        g = nx.convert_node_labels_to_integers(g)
        for key, val in all_bounds(g).items():
            if abs(val - 2 * k) > 1e-9:
                bad.append((key, f"{k}-regular", val, 2 * k))
        if bipartite and abs(laplacian_spectral_radius(g) - 2 * k) > 1e-9:
            bad.append(("mu", f"{k}-regular bipartite",
                        laplacian_spectral_radius(g), 2 * k))
    if verbose:
        print("selftest:", "PASS" if not bad else f"FAIL {bad}")
    return bad


if __name__ == "__main__":
    bad = selftest(verbose=True)
    print(f"{len(TYPE1) + len(TYPE2)} bounds transcribed, "
          f"{len(OPEN)} open, {len(REFUTED)} already refuted")
    raise SystemExit(1 if bad else 0)
