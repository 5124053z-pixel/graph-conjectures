"""The near-extremal family for bound 46, and its asymptotics.

The counterpart of `family44.py`. The quotient search kept returning 3-cell
matrices like

    [[0, 150, 0], [147, 0, 4], [0, 150, 0]]      degrees (150, 151, 150)

but cells 1 and 3 there carry identical parameters, so the shape collapses: it
is a **semiregular bipartite** graph with part degrees a and b. For such a graph
mu = a + b exactly, and the bound evaluates on the single edge type
(d_i, d_j, m_i, m_j) = (a, b, b, a), giving

    slack(a, b) = 2 + sqrt( 2(a^2 + b^2) - 16ab/(a+b) + 4 ) - (a + b).

At a = b it is exactly 0 (regular bipartite, where the whole family is tight).
The interesting direction is b = a - 1, and there the algebra is clean. Writing
a = k+1, b = k, the radicand is

    4k^2 - 4k + 2 + 4/(2k+1)  =  (2k-1)^2 + 1 + 4/(2k+1),

so sqrt(radicand) = (2k-1) + (1 + 4/(2k+1)) / (2(2k-1)) + O(k^-3) and

    slack  =  1/(4k)  +  O(k^-2),

positive for every k and tending to 0. Same conclusion as bound 44, at a much
slower rate: 44 approaches equality like a^-3, 46 like k^-1. That is exactly why
the searches bottomed out three orders of magnitude apart -- +1.05e-7 against
+2.02e-4 -- and it is a fact about the two bounds, not about the search.

Computed at 50 digits, which is more than the rate needs but keeps the same
discipline as family44.py.
"""
from __future__ import annotations

from mpmath import mp, mpf, sqrt

mp.dps = 50


def slack(a, b):
    """Slack of bound 46 on the semiregular bipartite graph with degrees a, b."""
    a, b = mpf(a), mpf(b)
    r = 2 * (a * a + b * b) - 16 * a * b / (a + b) + 4
    return 2 + sqrt(r) - (a + b) if r >= 0 else None


def main():
    print("bound 46 on semiregular bipartite graphs, degrees (a, b); mu = a + b\n")

    print("a = b (regular bipartite): the whole family is exactly tight")
    for k in (2, 10, 100):
        print(f"   a = b = {k:>4}   slack = {mp.nstr(slack(k, k), 8)}")

    print("\nb = a - 1, the direction the search actually descends:")
    print(f"{'k':>10}  {'slack':>22}  {'slack * 4k':>12}")
    for k in (10, 50, 150, 600, 5000, 100000, 10**7):
        s = slack(k + 1, k)
        print(f"{k:>10}  {mp.nstr(s, 12):>22}  {mp.nstr(s * 4 * k, 10):>12}")
    print("\n-> slack = 1/(4k) + O(k^-2); positive for every k, never crosses 0")

    print("\nwider gaps are worse, so b = a - 1 really is the extremal direction:")
    for gap in (1, 2, 3, 5, 10):
        s = slack(1000 + gap, 1000)
        print(f"   a - b = {gap:>2}   slack = {mp.nstr(s, 8)}")

    print("\ncompare bound 44 on its own family: slack = 1/(16 a^3).")
    print("46 approaches equality like k^-1, 44 like a^-3 -- which is why the")
    print("two searches bottomed out at +2.02e-4 and +1.05e-7 respectively.")


if __name__ == "__main__":
    main()
