"""Exhaustive generation of connected cubic graphs, count-verified.

McKay's public collections do not include a complete list of connected cubic
graphs, so they are generated here.

The rule: every connected cubic multigraph on n+2 vertices arises from one on n
vertices by choosing two edge slots, subdividing each, and joining the two new
vertices. (Choosing the same slot twice subdivides that edge twice and leaves
the two new vertices doubly joined.) Starting from the two connected cubic
multigraphs on two vertices -- three parallel edges, and two loops joined by an
edge -- this is applied repeatedly, and simple graphs are read off at each
level.

**Multigraph intermediates are not optional.** The first version of this file
used the same rule but kept only simple graphs, and produced 4 of the 5
connected cubic graphs on 8 vertices: one of them is reachable only through a
multigraph on 6 vertices. A silently incomplete enumeration turns "no
counterexample found" into a false statement, and is the one failure mode a scan
cannot detect from the inside -- so the generation rule is not trusted, it is
checked against the known counts (OEIS A002851) at every level, and the run
aborts on the first disagreement.

Note on isomorphism rejection: networkx's Weisfeiler-Lehman hash is useless
here. 1-WL never refines a regular graph, so every cubic graph on n vertices
gets the same hash. The adjacency spectrum is used instead, with an exact
isomorphism test inside each bucket.
"""
from __future__ import annotations

import argparse
import itertools
import time

import numpy as np
import networkx as nx

#: OEIS A002851, connected cubic *simple* graphs on n vertices.
KNOWN_COUNTS = {4: 1, 6: 2, 8: 5, 10: 19, 12: 85, 14: 509, 16: 4060,
                18: 41301, 20: 510489}

#: the two connected cubic multigraphs on two vertices
BASE = [
    [(0, 1), (0, 1), (0, 1)],          # theta
    [(0, 0), (1, 1), (0, 1)],          # dumbbell
]


def is_simple(edges):
    return (len(set(edges)) == len(edges)
            and all(u != v for u, v in edges))


def to_graph(edges, n):
    g = nx.Graph()
    g.add_nodes_from(range(n))
    g.add_edges_from(edges)
    return g


def spectrum_key(edges, n):
    """Isomorphism-invariant bucket key.

    Closed-walk counts rather than eigenvalues: `diag(A^k)` counts closed walks
    of length k at each vertex, is invariant under relabelling, and separates
    cubic graphs about as well as the spectrum at roughly a tenth of the cost.
    (`eigvalsh` was the whole bottleneck -- 80 microseconds per candidate, and
    there are millions of candidates.) Correctness does not depend on how good
    the key is: it only buckets, and an exact isomorphism test runs inside each
    bucket. A weak key costs time, never completeness."""
    a = np.zeros((n, n), dtype=np.int64)
    for u, v in edges:
        if u == v:
            a[u, u] += 2
        else:
            a[u, v] += 1
            a[v, u] += 1
    a2 = a @ a
    a3 = a2 @ a
    a4 = a2 @ a2
    return (tuple(sorted(np.diag(a2).tolist())),
            tuple(sorted(np.diag(a3).tolist())),
            tuple(sorted(np.diag(a4).tolist())),
            int(a4.sum()))


def multigraph(edges, n):
    g = nx.MultiGraph()
    g.add_nodes_from(range(n))
    g.add_edges_from(edges)
    return g


def connected(edges, n):
    seen, stack = {0}, [0]
    nbr = {v: [] for v in range(n)}
    for u, v in edges:
        nbr[u].append(v)
        nbr[v].append(u)
    while stack:
        for w in nbr[stack.pop()]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == n


def expand(edges, n):
    """All ways to grow a cubic multigraph by two vertices."""
    x, y = n, n + 1
    for i, j in itertools.combinations_with_replacement(range(len(edges)), 2):
        if i == j:
            a, b = edges[i]
            new = [e for k, e in enumerate(edges) if k != i]
            new += [(a, x), (x, y), (x, y), (y, b)]
        else:
            a, b = edges[i]
            c, d = edges[j]
            new = [e for k, e in enumerate(edges) if k not in (i, j)]
            new += [(a, x), (x, b), (c, y), (y, d), (x, y)]
        yield [tuple(sorted(e)) for e in new]


def generate_multigraphs(nmax, verbose=True):
    """{n: [edge lists of connected cubic multigraphs on n vertices]}.

    The intermediates of `generate`, exposed because the triangle-replacement
    family is indexed by them."""
    level = [list(map(lambda e: tuple(sorted(e)), b)) for b in BASE]
    out = {}
    for n in range(4, nmax + 1, 2):
        buckets = {}
        for edges in level:
            for cand in expand(edges, n - 2):
                if not connected(cand, n):
                    continue
                bucket = buckets.setdefault(spectrum_key(cand, n), [])
                gm = multigraph(cand, n)
                if not any(nx.is_isomorphic(gm, k[1]) for k in bucket):
                    bucket.append((cand, gm))
        level = [c for b in buckets.values() for c, _ in b]
        out[n] = level
        if verbose:
            print(f"  n = {n:>2}: {len(level):>9,} connected cubic multigraphs")
    return out


def generate(nmax, verbose=True, keep_multigraphs=True):
    """{n: [connected cubic simple graphs on n vertices]}."""
    level = [list(map(lambda e: tuple(sorted(e)), b)) for b in BASE]
    out = {}
    for n in range(4, nmax + 1, 2):
        t0 = time.time()
        buckets = {}
        for edges in level:
            for cand in expand(edges, n - 2):
                if not connected(cand, n):
                    continue
                key = spectrum_key(cand, n)
                bucket = buckets.setdefault(key, [])
                gm = multigraph(cand, n)
                if not any(nx.is_isomorphic(gm, k[1]) for k in bucket):
                    bucket.append((cand, gm))
        level = [c for b in buckets.values() for c, _ in b]
        simple = [to_graph(e, n) for e in level if is_simple(e)]
        out[n] = simple
        exp = KNOWN_COUNTS.get(n)
        ok = exp is None or len(simple) == exp
        if verbose:
            print(f"  n = {n:>2}: {len(simple):>8,} simple  "
                  f"(expected {exp if exp else '?':>8})   "
                  f"{len(level):>8,} multigraphs   "
                  f"{'OK' if ok else '*** MISMATCH ***'}   "
                  f"[{time.time() - t0:.1f}s]")
        if not ok:
            raise RuntimeError(
                f"generation is incomplete at n={n}: produced {len(simple)} "
                f"simple graphs, expected {exp}. The scan must not be trusted "
                f"until this is fixed.")
        if not keep_multigraphs:
            level = [e for e in level if is_simple(e)]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=14)
    args = ap.parse_args()
    print("generating connected cubic graphs, verifying counts against "
          "OEIS A002851\n")
    generate(args.n)
    print("\nall counts match; the enumeration is complete in this range")


if __name__ == "__main__":
    main()
