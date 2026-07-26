import networkx as nx, wowii_toolkit as T
from allgraphs import connected_graphs
WITNESS = {T.b_greedy_cds, T.b_greedy_path, T.b_exhibit_forest}
gs = connected_graphs(8)
print(f"\nstructural-bound coverage over {len(gs)} connected graphs on 2..8\n")
print(f"   {'conj':>5} {'residual n<=7':>14} {'residual n<=8':>14}")
for name, (bounds, rhs) in sorted(T.CONJ.items(), key=lambda t: int(t[0])):
    struct = [b for b in bounds if b not in WITNESS]
    r = sum(1 for g in gs if not any(b(g) >= rhs(g) for b in struct))
    print(f"   {name:>5} {'-':>14} {r:>14}", flush=True)
