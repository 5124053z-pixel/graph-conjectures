"""The near-extremal family for bound 44, and its asymptotics.

The 3-cell quotient search kept converging on one shape. Read off and cleaned
up, it is the one-parameter family

    B(a) = [[0, 0, a+3],
            [0, 0,   1],
            [a, 2,   0]]        degrees (a+3, 1, a+2)

-- a semiregular bipartite graph between cells 1 and 3, with two pendant
vertices hung on each cell-3 vertex. Along it the slack is positive but shrinks
like a^-3, so the family approaches equality without ever reaching it:

    slack(a)  =  1/(16 a^3)  +  O(a^-4).

This is computed at 60 decimal digits, because at a = 10^6 the slack is 6e-20
against quantities of size 10^6 -- some 26 orders of magnitude below double
precision, where a float64 answer would be pure noise. Compare failures 2 and 7
in the README: this family is exactly the regime where a naive computation
invents a counterexample.

The conclusion is negative but useful: the direction the search was heading is
closed, and refuting bound 44 requires leaving this family.
"""
from __future__ import annotations

from mpmath import mp, mpf, sqrt, polyroots

mp.dps = 60


def slack(a, b=2, c=None):
    """Exact-to-60-digits slack for B = [[0,0,c],[0,0,1],[a,b,0]]."""
    c = a + b + 1 if c is None else c
    a, b, c = mpf(a), mpf(b), mpf(c)
    s1, s2, s3 = c, mpf(1), a + b
    m1, m2 = s3, s3
    m3 = (a * s1 + b * s2) / s3

    M = [[s1, mpf(0), -c], [mpf(0), s2, mpf(-1)], [-a, -b, s3]]
    tr = M[0][0] + M[1][1] + M[2][2]
    m2c = ((M[0][0] * M[1][1] - M[0][1] * M[1][0])
           + (M[0][0] * M[2][2] - M[0][2] * M[2][0])
           + (M[1][1] * M[2][2] - M[1][2] * M[2][1]))
    det = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
           - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
           + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    mu = max(abs(x) for x in polyroots([mpf(1), -tr, m2c, -det],
                                       maxsteps=200, extraprec=300))

    def f(x, y, p, q):
        r = 2 * ((x - 1) ** 2 + (y - 1) ** 2 + p * q - x * y)
        return 2 + sqrt(r) if r >= 0 else mpf("-inf")

    return max(f(s1, s3, m1, m3), f(s2, s3, m2, m3)) - mu


def main():
    print("bound 44, family B(a) = [[0,0,a+3],[0,0,1],[a,2,0]]\n")
    print(f"{'a':>10}  {'slack':>26}  {'slack * a^3':>14}")
    for a in (10, 46, 100, 219, 500, 1000, 5000, 20000, 100000, 1000000):
        s = slack(a)
        print(f"{a:>10}  {mp.nstr(s, 12):>26}  {mp.nstr(s * a**3, 10):>14}")
    print("\n-> slack = 1/(16 a^3) + O(a^-4); positive for every a, never crosses 0")

    print("\nneighbouring parameters at a = 1000 (both b = 2 and c = a+b+1 "
          "are forced):")
    for b in (1, 2, 3):
        for delta in (0, 1, 2):
            s = slack(1000, b, 1000 + b + delta)
            print(f"  b={b}  c=a+b+{delta}   slack = {mp.nstr(s, 8)}")


if __name__ == "__main__":
    main()
