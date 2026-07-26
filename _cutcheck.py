import io, sys, math
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
import networkx as nx
from allgraphs import connected_graphs
import wowii as W

def avg_ecc(g):
    e = nx.eccentricity(g); return sum(e.values()) / g.number_of_nodes()

def maxcomp(g):
    return max(nx.number_connected_components(nx.restricted_view(g, [v], []))
               for v in g)

HYP = {
    '194':  lambda g: W.indep_number(g) <= 1 + W.avg_l(g),
    '198a': lambda g: W.largest_induced_bipartite(g) <= 2 + avg_ecc(g),
    '200':  lambda g: W.largest_induced_tree(g) == math.ceil(1 + W.avg_l(g)),
    '217':  lambda g: W.max_leaves_via_cds(g) <= 4 * (1 if W.residue(g) == 2 else 0) + 2,
}

print("claim to test: no graph with a vertex whose deletion leaves >= 3 components")
print("satisfies any of the four hypotheses. Such a graph has no Hamiltonian path,")
print("so one would be an immediate counterexample.\n", flush=True)

gs = connected_graphs(8, verbose=True)
gs = [g for g in gs if g.number_of_nodes() >= 3]
print(f"\n{len(gs):,} connected graphs on 3..8 vertices", flush=True)

bad = [g for g in gs if maxcomp(g) >= 3]
print(f"{len(bad):,} of them have a vertex whose deletion leaves >= 3 components\n", flush=True)

for nm, h in HYP.items():
    v = [g for g in bad if h(g)]
    tag = "OK -- claim holds" if not v else "*** COUNTEREXAMPLE ***"
    print(f"   #{nm:<5} satisfying the hypothesis: {len(v):>3}   {tag}", flush=True)
    for g in v[:2]:
        print(f"        {sorted(g.edges())}", flush=True)
