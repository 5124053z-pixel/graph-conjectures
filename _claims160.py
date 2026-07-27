"""Two claims that would finish the C4-free case of conjecture 160."""
import random
import networkx as nx
import wowii as W
from allgraphs import connected_graphs

ml = lambda g: max(W.local_independence(g, v) for v in g)
Tv = lambda g, v: g.subgraph(list(g[v])).number_of_edges()
mT = W.max_triangles_at_vertex
Ls = W.max_leaves_via_cds
D = lambda g: max(d for _, d in g.degree())
both = lambda g: any(W.local_independence(g, v) == ml(g) and Tv(g, v) == mT(g)
                     for v in g)

print("CLAIM A  C4-free  =>  max l + max T <= Delta + 1")
print("CLAIM B  C4-free, no vertex attains both maxima  =>  Ls >= Delta + 1\n")

gs = [g for g in connected_graphs(8) if not W.has_c4_subgraph(g)]
print(f"exhaustive, {len(gs)} C4-free graphs on 2..8:")
print(f"   A fails {sum(1 for g in gs if ml(g) + mT(g) > D(g) + 1)}")
print(f"   B fails {sum(1 for g in gs if not both(g) and Ls(g) < D(g) + 1)}",
      flush=True)

named = [("Petersen", nx.petersen_graph()), ("Heawood", nx.heawood_graph()),
         ("C_9", nx.cycle_graph(9)), ("P_12", nx.path_graph(12)),
         ("balanced tree 2,3", nx.balanced_tree(2, 3)),
         ("star K1,8", nx.star_graph(8)),
         ("friendship F3", nx.windmill_graph(3, 3)),
         ("triangle chain", nx.Graph([(0, 1), (1, 2), (0, 2), (2, 3), (3, 4),
                                      (2, 4), (4, 5), (5, 6), (4, 6)]))]
print(f"\n   {'graph':<18} {'maxl':>5} {'maxT':>5} {'Delta':>6} {'Ls':>4} "
      f"{'A':>5} {'both':>5} {'B':>5}")
for nm, g in named:
    g = nx.convert_node_labels_to_integers(g)
    if W.has_c4_subgraph(g):
        print(f"   {nm:<18}   (has a C4)")
        continue
    A = ml(g) + mT(g) <= D(g) + 1
    b = both(g)
    B = b or Ls(g) >= D(g) + 1
    print(f"   {nm:<18} {ml(g):>5} {mT(g):>5} {D(g):>6} {Ls(g):>4} "
          f"{str(A):>5} {str(b):>5} {str(B):>5}", flush=True)

print("\nrandom C4-free graphs, 9..13 vertices:")
rng = random.Random(3)
tot = fa = fb = 0
for n in range(9, 14):
    for _ in range(25):
        g = nx.gnp_random_graph(n, 0.25, seed=rng.randrange(10 ** 6))
        guard = 0
        while W.has_c4_subgraph(g) and g.number_of_edges() and guard < 200:
            g.remove_edge(*rng.choice(list(g.edges())))
            guard += 1
        if not nx.is_connected(g) or g.number_of_edges() == 0:
            continue
        tot += 1
        if ml(g) + mT(g) > D(g) + 1:
            fa += 1
        if not both(g) and Ls(g) < D(g) + 1:
            fb += 1
    print(f"   through n={n}: {tot} tested, A fails {fa}, B fails {fb}",
          flush=True)
