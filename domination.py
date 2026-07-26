"""Domination and total domination numbers, exact, with a validation suite.

Target conjecture (TxGraffiti; stated open in Davila's survey
[arXiv:2507.17780](https://arxiv.org/abs/2507.17780) and in
[arXiv:2409.19379](https://arxiv.org/abs/2409.19379)):

    if G and H are connected with n(G), n(H) >= 2, then
        gamma_t(G [] H)  >=  gamma(G x H)

where `[]` is the Cartesian product and `x` the direct (tensor) product.

Why this one. A search for work citing the conjecture turns up nothing at all --
no proofs, no partial results, no counterexample hunts. Compare the Laplacian
bounds, where two groups with deep-RL and Monte Carlo machinery had swept the
field and a third had closed it six weeks before this project started. Here the
field appears to be empty, and the natural first move -- checking every pair of
small graphs -- is one nobody seems to have made.

Both invariants are minimum set cover: gamma over the closed neighbourhoods
N[v], gamma_t over the open ones N(v). Products of two graphs on <= 7 vertices
have <= 49 vertices, which exact branch and bound handles comfortably.
"""
from __future__ import annotations

import math

import networkx as nx


def _cover_sets(g, closed):
    """Bitmask per vertex of the (closed or open) neighbourhood it covers."""
    idx = {v: i for i, v in enumerate(g)}
    sets = []
    for v in g:
        m = 1 << idx[v] if closed else 0
        for u in g[v]:
            m |= 1 << idx[u]
        sets.append(m)
    return sets, (1 << len(idx)) - 1


def min_set_cover(sets, universe):
    """Exact minimum number of sets covering `universe`. None if impossible."""
    sets = sorted(set(sets), key=lambda m: -bin(m).count("1"))
    if universe & ~_union(sets):
        return None
    biggest = bin(sets[0]).count("1")

    # greedy upper bound
    covered, ub = 0, 0
    while covered != universe:
        s = max(sets, key=lambda m: bin(m & ~covered).count("1"))
        if not (s & ~covered):
            return None
        covered |= s
        ub += 1
    best = ub

    # branch and bound: always branch on the hardest uncovered element
    covering = {}
    for i in range(universe.bit_length()):
        bit = 1 << i
        if universe & bit:
            covering[bit] = [s for s in sets if s & bit]

    def rec(covered, used):
        nonlocal best
        if covered == universe:
            best = min(best, used)
            return
        remaining = bin(universe & ~covered).count("1")
        if used + math.ceil(remaining / biggest) >= best:
            return
        target, options = None, None
        for bit, opts in covering.items():
            if covered & bit:
                continue
            live = [s for s in opts if s & ~covered]
            if options is None or len(live) < len(options):
                target, options = bit, live
                if len(live) <= 1:
                    break
        for s in sorted(options, key=lambda m: -bin(m & ~covered).count("1")):
            rec(covered | s, used + 1)

    rec(0, 0)
    return best


def _union(sets):
    u = 0
    for s in sets:
        u |= s
    return u


def domination_number(g):
    sets, universe = _cover_sets(g, closed=True)
    return min_set_cover(sets, universe)


def total_domination_number(g):
    """None when G has an isolated vertex, where gamma_t is undefined."""
    sets, universe = _cover_sets(g, closed=False)
    return min_set_cover(sets, universe)


def cartesian(g, h):
    return nx.convert_node_labels_to_integers(nx.cartesian_product(g, h))


def direct(g, h):
    return nx.convert_node_labels_to_integers(nx.tensor_product(g, h))


def slack(g, h):
    """gamma_t(G [] H) - gamma(G x H). Negative refutes the conjecture."""
    return total_domination_number(cartesian(g, h)) - domination_number(direct(g, h))


# ---------------------------------------------------------------------------
# validation
# ---------------------------------------------------------------------------

def _gamma_path_cycle(n):
    return math.ceil(n / 3)


def _gamma_t_path_cycle(n):
    return n // 2 + math.ceil(n / 4) - n // 4


#: The path and cycle rows are *computed* from the standard formulas
#: gamma(P_n) = gamma(C_n) = ceil(n/3) and
#: gamma_t(P_n) = gamma_t(C_n) = floor(n/2) + ceil(n/4) - floor(n/4),
#: not typed out. The first version typed them, got P_8 and C_8 wrong (5 instead
#: of 4, from mis-evaluating the formula by hand), and the suite reported two
#: MISMATCHES against correct code. A validation suite is code too, and
#: hand-transcribed expected values are exactly as untrustworthy as the
#: hand-transcribed formulas the rest of this repository refuses to trust.
CASES = [(f"P_{n}", nx.path_graph(n), _gamma_path_cycle(n), _gamma_t_path_cycle(n))
         for n in range(4, 10)]
CASES += [(f"C_{n}", nx.cycle_graph(n), _gamma_path_cycle(n), _gamma_t_path_cycle(n))
          for n in range(5, 10)]
#: these are literature values for named graphs, so they stay explicit
CASES += [
    ("K_5", nx.complete_graph(5), 1, 2),
    ("K_{3,3}", nx.complete_bipartite_graph(3, 3), 2, 2),
    ("K_{2,5}", nx.complete_bipartite_graph(2, 5), 2, 2),
    ("Petersen", nx.petersen_graph(), 3, 4),
    ("Q_3", nx.convert_node_labels_to_integers(nx.hypercube_graph(3)), 2, 4),
]


def validate(verbose=True):
    bad = []
    for name, g, gam, gam_t in CASES:
        g = nx.convert_node_labels_to_integers(g)
        d, t = domination_number(g), total_domination_number(g)
        ok = (d == gam) and (t == gam_t)
        if not ok:
            bad.append((name, d, gam, t, gam_t))
        if verbose:
            print(f"  {name:<10} gamma = {d} (want {gam})   "
                  f"gamma_t = {t} (want {gam_t})   {'OK' if ok else 'MISMATCH'}")
    return bad


if __name__ == "__main__":
    print("validating gamma and gamma_t against the standard formulas\n")
    bad = validate()
    print()
    if bad:
        print(f"MISMATCHES: {bad}")
        raise SystemExit(1)
    print("all values reproduced\n")
    print("the conjecture on a few pairs:")
    for a, b in [(nx.path_graph(3), nx.path_graph(3)),
                 (nx.cycle_graph(4), nx.cycle_graph(5)),
                 (nx.complete_graph(3), nx.path_graph(4)),
                 (nx.petersen_graph(), nx.path_graph(2))]:
        gt = total_domination_number(cartesian(a, b))
        gd = domination_number(direct(a, b))
        print(f"  |G|={len(a)} |H|={len(b)}:  gamma_t(G[]H) = {gt}   "
              f"gamma(GxH) = {gd}   slack = {gt - gd:+d}")
