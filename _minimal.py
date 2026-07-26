"""Smallest subset of structural bounds that settles each conjecture in range.

A pile of bounds covering a conjecture says little. The useful statement is
"these k inequalities imply it", which is a reduction -- as for conjecture 19,
where one bound sufficed and removed the bipartite number entirely.
"""
import itertools
import networkx as nx
import wowii as W
import wowii_toolkit as T
from allgraphs import connected_graphs

WIT = {T.b_greedy_cds, T.b_greedy_path, T.b_exhibit_forest}
NAME = {}
for n, b in [('acyclic', T.b_acyclic), ('star', T.b_star), ('diam+1', T.b_path),
             ('girth-1', T.b_cycle), ('caterpillar', T.b_caterpillar_cycle),
             ('cyclepath', T.b_cyclepath), ('2-step star', T.b_two_step_star),
             ('3-step star', T.b_three_step_star),
             ('1-vertex-off', T.b_one_vertex_off), ('2rad-1', T.b_rad_tree),
             ('cyclomatic', T.b_cyclomatic), ('cat-indep', T.b_caterpillar_indep),
             ('two-indep', T.b_two_indep), ('b/2+1', T.b_half_bipartite),
             ('rad path', T.b_rad_path), ('tree', W.largest_induced_tree),
             ('forest', W.largest_induced_forest)]:
    NAME[id(b)] = n

gs = connected_graphs(8)
N = len(gs)
print(f"{N} connected graphs on 2..8\n")

for name, (bounds, rhs) in sorted(T.CONJ.items(), key=lambda t: int(t[0])):
    st = [b for b in bounds if b not in WIT]
    r = [rhs(g) for g in gs]
    masks = []
    for b in st:
        m = 0
        for i, g in enumerate(gs):
            try:
                if b(g) >= r[i]:
                    m |= 1 << i
            except Exception:
                pass
        masks.append(m)
    full = (1 << N) - 1
    union = 0
    for m in masks:
        union |= m
    if union != full:
        miss = bin(full & ~union).count("1")
        print(f"   #{name:<4} structural residual {miss} -- no subset works",
              flush=True)
        continue
    best = None
    for k in (1, 2, 3, 4):
        for combo in itertools.combinations(range(len(st)), k):
            u = 0
            for i in combo:
                u |= masks[i]
            if u == full:
                best = combo
                break
        if best:
            break
    labs = ([NAME.get(id(st[i]), getattr(st[i], "__name__", "?")) for i in best]
            if best else None)
    print(f"   #{name:<4} needs {len(best) if best else '5+'}: {labs}", flush=True)
