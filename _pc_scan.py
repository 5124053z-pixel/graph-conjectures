import sys, random, time
import networkx as nx
from domination import total_domination_number as gt, domination_number as gd, cartesian, direct
from scan_products import connected_graphs
rng = random.Random(5)
Hs = [(f'P_{k}', nx.path_graph(k)) for k in range(3, 9)] + \
     [(f'C_{k}', nx.cycle_graph(k)) for k in range(3, 9)]
by = connected_graphs(6)
Gs = [g for n in sorted(by) for g in by[n] if n >= 4]
for n in range(7, 12):
    Gs += [nx.path_graph(n), nx.cycle_graph(n), nx.star_graph(n-1), nx.wheel_graph(n)]
    for _ in range(3):
        g = nx.gnp_random_graph(n, rng.uniform(0.2, 0.55), seed=rng.randint(0, 10**9))
        if nx.is_connected(g): Gs.append(g)
t0 = time.time()
for hn, h in Hs:
    worst, tight, tot = 9, 0, 0
    for g in Gs:
        if g.number_of_nodes() * h.number_of_nodes() > 40: continue
        s = gt(cartesian(g, h)) - gd(direct(g, h)); tot += 1
        if s == 0: tight += 1
        worst = min(worst, s)
        if s < 0:
            print(f'*** COUNTEREXAMPLE H={hn} |G|={g.number_of_nodes()} slack={s}', flush=True)
            print('   G:', sorted(g.edges()), flush=True)
    print(f'  H={hn:<4} {tot:>4} G tested  min slack {worst:+d}  '
          f'tight {tight:>4} ({100*tight/max(tot,1):>4.1f}%)  [{time.time()-t0:>4.0f}s]', flush=True)
