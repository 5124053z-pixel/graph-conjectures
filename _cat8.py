import io, sys, itertools
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import networkx as nx
from allgraphs import connected_graphs
import wowii as W

def girth(g): return 0 if nx.is_forest(g) else nx.girth(g)

def shortest_cycles(g):
    gi = girth(g)
    if gi == 0: return []
    out = []
    for s in itertools.combinations(list(g), gi):
        sub = g.subgraph(s)
        if sub.number_of_edges() == gi and all(d == 2 for _, d in sub.degree()) \
                and nx.is_connected(sub):
            out.append(set(s))
    return out

def caterpillar(g):
    gi = girth(g)
    if gi < 4: return 0
    best = 0
    for C in shortest_cycles(g):
        for x in C:
            P = C - {x}
            cand = [a for a in g if a not in C and len(set(g[a]) & P) == 1]
            m = W.indep_number(g.subgraph(cand)) if cand else 0
            best = max(best, (gi - 1) + m)
    return best

def ecc_c(g): return W.ecc_set(g, nx.center(g))
def r144(g): return girth(g) - 1 + ecc_c(g)

print("conjecture 144 over all connected graphs on <= 8 vertices,")
print("using only proved lower bounds\n", flush=True)
gs = connected_graphs(8, verbose=True)
gs = [g for g in gs if g.number_of_nodes() >= 3]
print(f"\n{len(gs):,} graphs", flush=True)

cheap = [lambda g: nx.diameter(g) + 1,
         lambda g: 1 + max(W.local_independence(g, v) for v in g),
         lambda g: max(girth(g) - 1, 0)]

stage1 = [g for g in gs if not any(f(g) >= r144(g) for f in cheap)]
print(f"not settled by the three cheap bounds : {len(stage1)}", flush=True)

stage2 = [g for g in stage1 if caterpillar(g) < r144(g)]
print(f"still not settled by the caterpillar   : {len(stage2)}", flush=True)

bad = [g for g in gs if W.largest_induced_tree(g) < r144(g)]
print(f"actual violations of conjecture 144    : {len(bad)}", flush=True)
for g in stage2[:5]:
    print(f"   residual: n={g.number_of_nodes()} girth={girth(g)} "
          f"rhs={r144(g)} tree={W.largest_induced_tree(g)} "
          f"cat={caterpillar(g)}  {sorted(g.edges())}", flush=True)
