# Counterexample search for graph-invariant inequalities

Status: amateur/independent investigation. Started 2026-07-26.

## Why this problem class

Most conjectures of this kind have a single shape:

```
for every (connected) graph G on n vertices,   f(G) ≥ g(G)
```

with `f, g` computable graph invariants. That makes the whole workflow mechanical:

* a **counterexample is one finite graph**, small enough to write down;
* **verification is a calculation**, with no human judgement involved;
* the **search is "minimise the slack"** `f(G) − g(G)`, so switching conjectures
  means changing one function and nothing else.

This is the structure exploited in Adam Zsolt Wagner, *Constructions in
combinatorics via neural networks* ([arXiv:2104.14516](https://arxiv.org/abs/2104.14516)),
which refuted several open conjectures this way, and in the follow-up
literature (Ghebleh et al. 2024; *Parallelizing Wagner's approach*,
[arXiv:2509.01607](https://arxiv.org/pdf/2509.01607)).

Material is not scarce: conjecture-generating software such as **AutoGraphiX**
and Fajtlowicz's **Graffiti** has produced large numbers of candidate
inequalities, and the backlog of unresolved ones is the raw material here.

## Method note: validate the instrument before pointing it at an open problem

The first thing built here was run against a conjecture whose answer is
**already known**, so that a wrong tool would show up as a disagreement with the
literature rather than as a fake discovery. This is deliberate; see "What went
wrong" below for two ways the naive approach failed.

## Validation target

**Conjecture** (Aouchiche et al., originally produced by AutoGraphiX; Conjecture 2.1
in Wagner's paper). *For every connected graph `G` on `n ≥ 3` vertices with largest
adjacency eigenvalue `λ₁` and matching number `μ`:*

```
λ₁ + μ  ≥  √(n−1) + 1
```

Known to be false: Stevanović gave a 600-vertex counterexample; Wagner's search
found one on 19 vertices.

### Result: reproduced, and the minimum is 18 vertices, not 19

Exhaustive enumeration of **all non-isomorphic trees** settles this completely
(`exhaustive_trees.py`):

| n | trees | min slack | counterexamples |
|---|---|---|---|
| ≤13 | | `+0.000000000000` | 0 |
| 14 | 3,159 | `−8.9e-16` | 0 — **spurious**, see below |
| 15–17 | | `+0.000000000000` | 0 |
| **18** | 123,867 | **`−0.021810091692`** | **1** |
| 19 | 317,955 | `−0.080363026951` | 2 |

The 18-vertex counterexample has degree sequence `[9, 8, 2, 1×15]` — two star
centres joined through a single degree-2 vertex — with

```
λ₁ = 3.101295…   (λ₁² = 8 + φ, φ the golden ratio)
μ  = 2
λ₁ + μ = 5.101295… < √17 + 1 = 5.123105…
```

**Restricting to trees loses nothing.** Adding an edge to a connected graph
strictly increases `λ₁` (Perron–Frobenius) and cannot decrease `μ`, so the
minimum of `λ₁ + μ` over connected graphs on `n` vertices is always attained on a
spanning tree. The table above is therefore a definitive statement about *all*
connected graphs, not just trees: **no counterexample exists below n = 18.**

Wagner reports a 19-vertex example but makes no minimality claim, so this is a
small refinement rather than a correction.

## What went wrong (recorded rather than hidden)

**1. The stochastic search failed; brute force succeeded.** `search.py` implements
the cross-entropy method over the `n(n−1)/2` edge indicators — the canonical
Wagner setup. At `n = 19` it plateaued at slack `+3.48` and never came close.
The reason is visible in hindsight: the star `K₁,ₙ₋₁` attains the bound with
*equality*, so counterexamples live in a narrow region near a very sparse graph,
while the search starts at edge-probability 0.5 (average degree ≈ 9). Exhaustive
enumeration of all 317,955 trees found both counterexamples in 38 seconds.
For `n` in this range, brute force is simply better, and the fancy method was
the wrong instrument.

**2. A floating-point false positive at n = 14.** The first scan reported a
counterexample at `n = 14`. It is the star `K₁,₁₃`, for which `λ₁ = √13` and
`μ = 1`, so `λ₁ + μ` equals the threshold *exactly*; the computed `λ₁²` came out
as `12.999999999999996` and the slack as `−8.9 × 10⁻¹⁶`. Reporting `n = 14` as
the minimum would have been a fabricated result caused entirely by rounding.
All counts above therefore use a tolerance of `slack < −10⁻⁹`, and the spurious
case is listed separately rather than dropped silently.

The general lesson, which applies to the whole problem class: **these conjectures
are frequently tight**, with the extremal graph achieving equality. Any search
that reports a counterexample with slack at the 10⁻¹⁵ level has found rounding
error, not mathematics. Exact or interval arithmetic is needed near equality.

---

# Open problem: 36 conjectured bounds on the Laplacian spectral radius

Started 2026-07-26, immediately after the validation above.

## The target

Let `μ(G)` be the largest eigenvalue of the Laplacian `D − A`, let `d_v` be the
degree of `v`, and let `m_v` be the average degree of the neighbours of `v`.
There is a family of 68 automatically generated conjectures of the form

```
μ(G) ≤ max_v f(d_v, m_v)              ("type 1")
μ(G) ≤ max_{uv ∈ E} f(d_u, d_v, m_u, m_v)   ("type 2")
```

They were produced by Brankov, Hansen and Stevanović from a database of 273,214
connected graphs on at most 9 vertices, and they have since been attacked twice:

| attack | coverage | result |
|---|---|---|
| Al-Yakoob–Ghebleh–Kanso–Stevanović, deep cross-entropy ([arXiv:2509.01607](https://arxiv.org/abs/2509.01607)) | fixed sizes n = 20–24; plus exhaustive over subquartic graphs, n ≤ 14 | 30 refuted |
| Taieb–Roucairol–Cazenave–Harutyunyan, 8 Monte Carlo search algorithms ([LION 2025](https://www.lamsade.dauphine.fr/~cazenave/papers/ArticleLioraLNCS.pdf)) | incremental construction to n ≈ 20–30 | 29 refuted, of which 2 (45, 48) were previously open |

**36 remain open.** They satisfy the selection criteria above: `μ`, `d_v`, `m_v`
are all polynomial-time, and switching between conjectures is a one-line change.

`laplacian_bounds.py` transcribes all 68 from the LION appendix.

### Guarding the transcription

The formulas were read off a PDF, so they are checked rather than trusted. All
68 bounds share an identity: **on a `k`-regular graph every one of them
evaluates to exactly `2k`**, and on a `k`-regular *bipartite* graph `μ = 2k`
exactly. Any transcription slip breaks that identity. `selftest()` runs the
check on `K_{3,3}`, `K_{4,4}`, `C₈`, `Q₃` and the Petersen graph before any
scan, and all 68 pass.

This is the preventive form of the `n = 14` floating-point lesson from the
validation phase: the whole family is *exactly tight* on regular bipartite
graphs, which is precisely the regime where a small error is invisible.

### Conjectures 10 and 23 are the same conjecture

Bound 10 is `√(d(d + 3m))` and bound 23 is `√(d² + 3dm)`. These are the same
expression. Both are listed separately, and both are listed as open. Checked by
evaluating all 68 formulas at 4,000 random points and comparing (`duplicates`
check); it is the only coincidence among the 68. So the number of *distinct*
open conjectures is 35, not 36.

## Result: 12 of the 36 are true, and provable in one line each

Everyone attacking these conjectures has been *searching for counterexamples* —
deep reinforcement learning, Monte Carlo tree search, exhaustive generation. The
opposite direction appears not to have been tried: checking whether a
conjectured bound is simply **implied by a theorem that has been known since the
1980s**.

Because all these expressions are local, the implication reduces to a pointwise
inequality between two small algebraic functions. `dominance.py` scans for
violations of that inequality over the feasible region; where none is found, the
implication is then verified by hand. Seven survived:

| conj | conjectured bound `μ ≤ …` | implied by | pointwise inequality |
|---|---|---|---|
| 5 | `max_v (d²/m + m)` | Anderson–Morley | `d²/m + m − 2d = (d−m)²/m ≥ 0` |
| 6 | `max_v √(m² + 3d²)` | Li–Zhang | `(m²+3d²) − 2d(d+m) = (d−m)² ≥ 0` |
| 9 | `max_v (m + 3d)/2` | Li–Zhang | `((m+3d)/2)² − 2d(d+m) = (d−m)²/4 ≥ 0` |
| 12 | `max_v √(2m² + 2d²)` | Merris | `2m²+2d² − (d+m)² = (d−m)² ≥ 0` |
| 34 | `max_{uv∈E} 2(d_u²+d_v²)/(d_u+d_v)` | Anderson–Morley | `2(a²+b²) − (a+b)² = (a−b)² ≥ 0` |
| 37 | `max_{uv∈E} √(2(d_u²+d_v²))` | Anderson–Morley | same |
| 38 | `max_{uv∈E} [2 + √(2(d_u−1)² + 2(d_v−1)²)]` | Anderson–Morley | `√(2(x²+y²)) ≥ x+y`, `x=d_u−1` |

Every one of them is the *same* inequality — `(a−b)² ≥ 0`, i.e. QM ≥ AM — wearing
a different disguise.

### Five more, by pinning the vertex instead of quantifying over all of them

Pointwise domination asks for `B_v ≥ M_v` at **every** vertex, which is far more
than is needed: only `max_v B_v ≥ μ` is required. The gap closes by naming a
vertex where the inequality is easy.

> Let `v*` have maximum degree `Δ`. Every neighbour of `v*` has degree `≤ Δ`, and
> `m_{v*}` is an average of those degrees, so **`m_{v*} ≤ Δ = d_{v*}`: at a vertex
> of maximum degree, `d ≥ m` always.** Anderson–Morley gives `μ ≤ 2Δ`. Hence if
> `B(d, m) ≥ 2d` merely on the half-plane `m ≤ d`, then
> `max_v B_v ≥ B_{v*} ≥ 2Δ ≥ μ`.

That condition is much weaker than pointwise domination, and five more
conjectures satisfy it — all of them collapsing to `d ≥ m` itself:

| conj | bound | reduces to |
|---|---|---|
| 1 | `max_v √(4d³/m)` | `4d³/m ≥ 4d²` ⟺ `d ≥ m` |
| 4 | `max_v 2d²/m` | `2d²/m ≥ 2d` ⟺ `d ≥ m` |
| 7 | `max_v (d²/m + d)` | `d²/m ≥ d` ⟺ `d ≥ m` |
| 14 | `max_v 2d³/m²` | `d³/m² ≥ d` ⟺ `d ≥ m` |
| 16 | `max_v 2d⁴/m³` | `d⁴/m³ ≥ d` ⟺ `d ≥ m` |

Conjecture 5 also falls to this argument independently: `d²/m + m ≥ 2d` by AM–GM
holds at every vertex, so the bound is `≥ 2Δ` outright.

Theorems used, all classical:

* **Anderson–Morley** (1985): `μ ≤ max_{uv∈E} (d_u + d_v)`
* **Merris** (1998): `μ ≤ max_v (d_v + m_v)`
* **Li–Zhang** (1998): `μ ≤ max_v √(2 d_v (d_v + m_v))`
* **Li–Pan** (2000): `μ ≤ max_v (d_v + √(d_v² + 8 d_v m_v))/2`
* **Das** (2003): `μ ≤ max_{uv∈E} (d_u(d_u+m_u) + d_v(d_v+m_v))/(d_u+d_v)`

`known_bounds.py` re-validates all five numerically on 7,500 graphs before they
are used. The corpus is used as a **falsifier, not a certifier** — it can expose
a misremembered theorem, but it cannot license one, so only bounds with a
citation are admitted. `proofs.py` re-checks each of the seven pointwise
inequalities as a regression test.

**Status: 36 open → 24 open** (23 distinct, since 10 and 23 are the same
conjecture and both remain).

Resolved here: **1, 4, 5, 6, 7, 9, 12, 14, 16, 34, 37, 38**.

### …all twelve of which had been published six weeks earlier

**Damnjanović, Ha and Stevanović, *Upper bounds for the Laplacian spectral
radius: proofs and counterexamples*, [arXiv:2606.14550](https://arxiv.org/abs/2606.14550),
12 June 2026.** Of the 36, they confirm 22, refute 12, and leave **two** open.

Every one of the twelve above is in their confirmed list, and the proofs are not
merely equivalent but identical:

* their Proposition 3.1 is the max-degree argument (`m_i ≤ Δ = d_i`, `μ ≤ 2Δ`)
  applied to bounds 1, 4, 5, 7, 14, 16 — the same six, by the same route;
* their Proposition 3.2 handles 6 and 9 via Li–Pan, 3.3 handles 12 via Merris,
  3.4 handles 34, 37, 38 via Anderson–Morley.

They also record that bounds 10 and 23 are identical.

So the correct account of this section is: **nothing here is new.** What it is
instead is an unusually strong validation of the method — an independent
re-derivation, from a PDF table of formulas, of a published classification,
agreeing on all twelve results and on the duplicate. The tooling works. It was
pointed at ground that had already been cleared.

### The two that are actually open

| conj | bound |
|---|---|
| **44** | `μ ≤ max_{ij∈E} [ 2 + √(2((d_i−1)² + (d_j−1)² + m_i m_j − d_i d_j)) ]` |
| **46** | `μ ≤ max_{ij∈E} [ 2 + √(2(d_i² + d_j²) − 16 d_i d_j/(m_i + m_j) + 4) ]` |

These are the only candidate bounds in the family whose status is unknown, and
the state of the art on them is six weeks old. Everything below targets them.

Convention (theirs, adopted here): if the radicand is negative the term is `−∞`,
i.e. that edge is dropped from the maximum rather than clamped. This makes the
bound *smaller* and so is the refuter-friendly reading; clamping to 0 instead,
as the first implementation here did, silently makes the conjecture harder to
refute.

## Attacking 44 and 46: search over equitable partitions

`quotient_search.py`. The method is the one that refuted the other twelve, and it
rests on a simplification worth stating plainly.

If `G` has an equitable partition into cells `C_1…C_k` with quotient matrix
`B = (b_ij)` — every vertex of `C_i` has exactly `b_ij` neighbours in `C_j` — then
every vertex of `C_i` has the same degree `s_i = Σ_j b_ij` and the same
neighbour-average `m_i = (Σ_j b_ij s_j)/s_i`. So the conjecture's right-hand side
is **a function of `B` alone**. And `μ(G) ≥ ρ(L_B)` where `L_B = diag(s) − B`.
Hence

```
ρ(L_B)  >  max_{b_ij > 0} f(s_i, s_j, m_i, m_j)
```

certifies a counterexample — with no graph ever constructed and no bound on how
large it has to be. Cell sizes enter only through realizability, and since
multiplying every `n_i` by a common even factor preserves `n_i b_ij = n_j b_ji`
while making the other conditions easier, **realizability reduces to the ratio
system being consistent**. Parametrising by cell sizes makes that automatic:
`b_ij = c·n_j/gcd(n_i,n_j)`, `b_ji = c·n_i/gcd(n_i,n_j)`.

So the search space is small integer matrices, and it covers graphs of unbounded
size. Simulated annealing over `(n, c)`.

### The instrument is validated against twelve known answers

Table 3 of arXiv:2606.14550 gives the quotient matrices refuting bounds 11, 13,
18–22, 24, 30, 40, 47 and 56. `validate()` recomputes all twelve and requires
each to come out as a counterexample; this exercises the quotient arithmetic,
the `−∞` convention and the bound transcriptions simultaneously. All twelve
reproduce, with slacks from `−0.0057` to `−1.27`.

That is a far stronger check than any self-test: twelve independent known
answers, none of which the code was written to produce.

### Status: no counterexample, and one direction closed

Nothing refuted. What has been established is where the counterexamples are
*not*.

**Exhaustive over all graphs on ≤ 7 vertices and all trees on ≤ 18 vertices**
(205,997 graphs): no counterexample to either bound.

**Exhaustive over all 2-cell quotients** with diagonal entries `< 60` and
off-diagonal `< 120`: minimum slack exactly `0` for both bounds, attained at
`B = [[0,1],[1,0]]` — a single edge, i.e. the regular bipartite case again.
Since a 2-cell quotient has `L_B = [[b,−b],[−c,c]]` and therefore `μ ≥ b + c`
exactly, this family is small enough to settle: **no 2-cell counterexample
exists to either bound.** Three or more cells are necessary.

### Bound 44 is asymptotically tight along an explicit family

The 3-cell search kept converging on one shape, and it turned out to be a
one-parameter family. With

```
B = [[0, 0, a+3],
     [0, 0,   1],
     [a, 2,   0]]        degrees (a+3, 1, a+2)
```

— a semiregular bipartite graph between cells 1 and 3, with two pendant vertices
per cell-3 vertex — the slack is positive but shrinks cubically. Computed at 60
decimal digits:

| `a` | slack | `slack · a³` |
|---|---|---|
| 1,000 | `6.2313e-11` | 0.0623129 |
| 5,000 | `4.9970e-13` | 0.0624625 |
| 20,000 | `7.8113e-15` | 0.0624906 |
| 100,000 | `6.2498e-17` | 0.0624981 |
| 1,000,000 | `6.2500e-20` | 0.0624998 |

So along this family

```
slack  =  1/(16 a³)  +  O(a⁻⁴)
```

positive for every `a`, tending to 0. Neighbouring parameters are worse by
orders of magnitude — `b = 2` and `c = a + b + 1` are both forced; changing `c`
by one costs a factor of 10⁷.

This does not prove bound 44, but it closes the direction the search was
heading: the most promising numerical family is provably (to 60 digits, out to
`a = 10⁶`) on the wrong side, with a clean asymptotic that never crosses zero.
Anyone refuting bound 44 has to leave this family.

### Bound 46 is asymptotically tight too, along a different family and at a different rate

The same treatment applied to bound 46. The quotient search kept returning
3-cell matrices like `[[0,150,0],[147,0,4],[0,150,0]]` with degrees
`(150, 151, 150)` — but cells 1 and 3 there carry identical parameters, so the
shape collapses to a **semiregular bipartite** graph with part degrees `a` and
`b`. There `μ = a + b` exactly and the bound has a single edge type, giving

```
slack(a,b) = 2 + √( 2(a²+b²) − 16ab/(a+b) + 4 ) − (a+b).
```

At `a = b` this is exactly 0 — regular bipartite, where the whole family is
tight. The descent direction is `b = a − 1`, and there the radicand is exactly
`(2k−1)² + 1 + 4/(2k+1)` for `a = k+1`, `b = k`, giving

```
slack = 1/(4k) + O(k⁻²),
```

positive for every `k`, tending to 0. Confirmed at 50 digits out to `k = 10⁷`,
where `slack · 4k = 1.00000025`. Widening the degree gap makes it worse, so
`b = a − 1` is the extremal direction.

**This also explains the searches.** Bound 44 approaches equality like `a⁻³` and
bound 46 like `k⁻¹`, so the two annealing runs bottoming out three orders of
magnitude apart — `+1.05e-7` against `+2.02e-4` — is a fact about the bounds,
not about the search. The best quotients found match the family values exactly:
`k = 150` gives `1.694e-3` and the annealer found `1.694e-3`; `k = 600` gives
`4.18e-4` against `4.19e-4`.

Both remaining bounds therefore have an explicit family approaching equality
from above and never crossing. Refuting either means leaving the family the
search descends into.

### Bound 46: the edge-dropout mechanism is a dead end

Bound 46 has a feature none of the others do — its radicand

```
R(i,j) = 2(d_i² + d_j²) − 16 d_i d_j/(m_i + m_j) + 4
```

can go **negative**, and under the paper's convention such an edge is dropped
from the maximum entirely. Dropping the edge that would have been the argmax
shrinks the bound, which looks like a purpose-built refutation route. It is not.

For `d_i = d_j = k`, `R < 0` iff `m_i + m_j < 4k²/(k²+1)`. Since `i` and `j` are
adjacent, `m_i ≥ (d_j + d_i − 1)/d_i = 2 − 1/k`, so `m_i + m_j ≥ 4 − 2/k`, and
`4 − 2/k < 4k²/(k²+1) ⟺ (k−1)² > 0`. A window therefore exists for every
`k ≥ 2`, and it is attained exactly at the **balanced double star** — two
adjacent centres of degree `k`, every other neighbour a leaf — where
`m_i = m_j = 2 − 1/k` and

```
R = −4(k−1)²/(2k−1)  <  0.
```

But forcing `m_i` that low *requires* the leaves, and a leaf edge never drops:
for a leaf `w` on a vertex `u` of degree `d ≥ 2`,
`R(w,u) = 2d² + 6 − 16d/(d + m_u) > 2d² − 10 > 0`. The surviving leaf edge
contributes `≈ 2 + d√2`, while `μ` of a double star is `≈ d + 1`. Dropping the
high-degree edge buys nothing, because the structure that makes it drop
introduces edges with larger values still.

Empirically: of the 205,997 graphs scanned, **1,496 have at least one dropped
edge and none has all edges dropped**, and none comes near refuting.

## Trying to prove 44 and 46 instead of refuting them

After several hours of search with nothing below `+1.05e-7` (bound 44) and
`+2.5e-4` (bound 46), and with bound 44's most promising family provably
positive, the sensible move is to try to prove them. Three attempts, all
recorded because all failed.

**The Perron bound narrows the gap but does not close it.** The proof of
Proposition 3.5 in arXiv:2606.14550 contains, as a step rather than a statement,
the bound

```
μ(G) ≤ max_{ij∈E} ( d_i + d_j + √((d_i−d_j)² + 4 m_i m_j) ) / 2
```

obtained from a Perron vector of `P⁻¹QP`. Added to `known_bounds.py` (and
re-validated on 7,500 graphs), it is by far the strongest of the six for this
purpose: it cuts bound 44's pointwise deficit from `−11.74` (the next best, Das)
to **`−0.79`**. Still negative, so no pointwise domination.

**Where it fails says why a pointwise proof cannot work.** The violations sit
exactly at the *minimum feasible* `m`: for `d_i = d_j = k` the worst point is
`m_i = m_j = 2 − 1/k`, which is the double-star configuration where every other
neighbour is a leaf. But bound 44 is a maximum over edges, and in a double star
the *leaf* edges are large — at `k = 2` the leaf edge gives 4 against the Perron
value 3.5. The bound survives through an edge other than the one where the
pointwise test is applied, so any proof has to be global.

**Lemma 3.6 does not reach edge-indexed bounds.** `collatz_wielandt.py` rebuilds
the machine behind most of the paper's confirmations: for `λ > Δ` and concave
`φ`, if `d_i φ(m_i/λ) ≤ (λ − d_i) φ(d_i/λ)` at every vertex then `μ ≤ λ`. For a
power `φ(x) = x^α` and a type-1 bound `B(d,m) = d·f(r)` with `r = m/d`, this
collapses to the one-variable condition `f(r) ≥ 1 + r^α`.

The machine is validated against the paper before being pointed anywhere: it
reproduces their `α = 1/4` for bound 8 with margin exactly 0, it fails on bound
27 at exactly the place where their own proof needs an extra argument, and it
finds **no exponent at all for any of the nine refuted type-1 bounds** — the
cheapest possible check that it cannot manufacture a proof.

Extending it to 44 and 46 needs a lower bound on `λ = max_{ij∈E} B(i,j)` in terms
of `(d_i, m_i)` alone. The only thing guaranteed at a vertex is that *some*
neighbour has degree `≥ m_i`, and that neighbour's own `m_j` is at least
`(d_i + d_j − 1)/d_j`. That is far too weak: the resulting condition fails by
`−14` to `−17` for bound 44 and `−2` to `−11` for bound 46, across every
exponent tried. The guarantee about one neighbour says almost nothing about a
maximum over all edges.

Every remaining open bound in the family is edge-indexed, and this is presumably
why they are the two that are left.

## Status of all 22, after the free-bound pass

`wowii_toolkit.py` collects every bound and reports the residual each conjecture
is left with over the 995 connected graphs on `≤ 7` vertices.

**Proved outright, for every graph:**

| conj | how |
|---|---|
| **65** | `distMin` only takes `0` and `1` in a connected graph, so the RHS never exceeds 2 |
| **141** | neighbourhood star for `girth ≤ 5`, locally-tree-like ball for `girth ≥ 6` |
| **316** | the hypothesis admits only complete graphs and `K₁,K₂,K₃` with pendants |
| **322** | the hypothesis forces `G = K_n` |

**Settled over every connected graph on `≤ 7` vertices by proved lower bounds:**

| conj | statement |
|---|---|
| **19** | `⌊avg ecc + max_v l(v)⌋ ≤ b(G)` — closed by Lemma B, from a residual of 148 |
| **40** | `f(G) ≥ ⌈(p+b+1)/2⌉` — closed by Lemma C; also proved outright for every graph with a Hamiltonian path |
| **61** | `f(G) ≥ residue + ⌈diam/3⌉` — closed by Lemma C; also proved outright for `diam ≤ 3` |
| **144** | `girth − 1 + ecc(centre) ≤ tree(G)` — closed by the cycle-path and caterpillar bounds together |

**Small residuals left:**

| conj | residual | note |
|---|---|---|
| 133 | **1** | proved outright for every graph containing a 4-cycle (905 of 995) |
| 142 | **2** | |
| 59 | **2** | down from 22; replacing `residue` by `α` makes it **false**, so `residue` is essential |
| 146 | 8 | all eight tight, all with `tree = n − 1` |

Four more bounds found in the last pass, all with zero violations:
`path(G) ≥ 2·rad − 1`; `tree(G) ≥ 2 + max_v l(v)` when `girth ≥ 5`;
`f(G) ≥ n − (m − n + 1)`, deleting one vertex per independent cycle; and the
containments `f ≥ tree ≥ path`. The cyclomatic one alone took `#59` from 22 to 2.

**A note on consolidating.** When the toolkit was tidied into `wowii_toolkit.py`
the single-path cycle bound was dropped in favour of the caterpillar version,
on the assumption that the more general one subsumed it. It does not: the
caterpillar attaches **single vertices at several** cycle vertices, the other
attaches **a whole path at one**, and neither contains the other. Conjecture 144
regained a residual until both were present — its one hard case is `C₅` with a
two-vertex path hanging off a single cycle vertex, which only the path version
sees.

**Elsewhere:** `160` is false as formalised (five-vertex counterexample); `194`,
`198a`, `200` and `217` are 80–99% covered by Chvátal–Erdős, with `217` split so
that the `residue ≠ 2` half is trivial; `314` sits 35/40 inside a published
characterisation; `2`, `100`, `145`, `291` remain open with no reduction.

None of the bounds is deep. What they are is **unused** — these statements had
been on DeLaVina's list since Graffiti.pc produced them, and nobody had checked
which of them a star, a shortest path or a chordless cycle already settles.

## What went wrong (continued)

**3. The tightness ranking was worthless.** The plan was to rank the 36 open
conjectures by minimum slack over small graphs, on the reasoning — correct in
the validation phase — that tight conjectures are the likely-false ones. A scan
over all connected graphs on ≤ 7 vertices and all trees on ≤ 18 vertices
returned *slack exactly 0 for all 36*. The reason is structural and should have
been predicted: every bound equals `2k` on a `k`-regular graph and `μ = 2k` on a
regular bipartite one, so `C₄` alone makes all 36 tight. Tightness carries no
information in this family, and the whole ranking idea was discarded.

The recoverable lesson: "tight conjectures are likely false" needs the tight
case to be *isolated*. Here equality holds on an infinite family that is
plainly not extremal for anything, so the heuristic degenerates.

**4. The cheap direction was skipped by everyone — including, at first, by me.**
These conjectures had resisted deep reinforcement learning and eight Monte Carlo
search algorithms, not because they are hard, but because those methods can only
ever answer "no counterexample found". Checking a conjectured bound against
theorems from the 1980s costs nothing and settles twelve of them.

**5. The literature check was one paper out of date, and that cost the whole
result.** Before starting, the status of these conjectures was established from
Ghebleh et al. (2026) and Taieb et al. (LION 2025), both of which say 36 remain
open. That was taken as current. It was six weeks stale:
[arXiv:2606.14550](https://arxiv.org/abs/2606.14550) had already reduced 36 to
2, by the same arguments, in June 2026.

The specific mistake is narrow and worth naming precisely. The searches that
were run were for *the conjectures* and for *the attacks on them* — Graffiti,
AutoGraphiX, Wagner, the Monte Carlo papers. No search was run for **recent work
citing the conjecture list itself**, which is where a resolution would appear.
An arXiv listing search on the source paper's citations would have surfaced it
in one query, before any code was written.

The general form: in a field where the answer may have arrived last month, the
question "is this open?" has to be asked *of the most recent citing work*, not
of the most recent paper one happens to have found. Checking the status of an
open problem is itself a step that can be done wrong.

**6. The annealer spent two runs stuck at slack exactly 0, twice for different
reasons.** First it collapsed onto *regular* quotients, where every bound in the
family equals `2s` and `μ ≤ 2s`, so slack is identically 0 — an enormous flat
plateau with no gradient out. Excluding regular quotients did not help: it then
found *disconnected* quotients containing one regular bipartite component, which
pins the maximum at exactly 0 regardless of the rest of the matrix.

Both are the same failure as #3 in a new costume — this family is exactly tight
on regular bipartite graphs, and any search that can reach them will sit there.
The fix is to prune both cases, which loses nothing (neither can be a
counterexample, and any counterexample can be taken connected). Only after both
prunings did the search start returning informative numbers at all.

The general lesson: when a conjecture family has a large exact-equality set, the
equality set is not a promising region to search — it is a trap that absorbs the
optimiser. It has to be excluded explicitly before search is meaningful.

**7. The exhaustive scan reported K₃,₃ as a counterexample, which is failure 2
verbatim.** The scan over 205,997 graphs used `slack < 0` as the test and duly
printed a "counterexample" to both bounds: the 3-regular bipartite graph K₃,₃,
where `μ = 6` exactly and the bound is exactly `6`, computed as `−1.1e-16`.

This is the *same* mistake as the `n = 14` star in the validation phase — the
one already written down at the top of this README as the lesson of the whole
project. Having documented it did not prevent repeating it in a new file three
commits later. The tolerance is now a shared constant rather than a literal
retyped per scan, which is the only version of the fix that actually holds.

**8. The one bound-46-specific idea did not survive contact with the algebra.**
Bound 46's radicand can go negative, dropping edges out of the maximum — a
mechanism no other bound in the family has, and the obvious place to look. Two
lines of algebra (above) show the structure required to make a high-degree edge
drop necessarily introduces leaf edges whose surviving value is larger than what
was dropped. Worth recording because it looked, for about an hour, like the
whole answer.

---

# Second open problem: zero forcing vs independence in cubic graphs

Started 2026-07-26, after the inventory below.

## The target

**Conjecture** (TxGraffiti, open since 2017; Davila et al.,
[arXiv:2410.21724](https://arxiv.org/abs/2410.21724)). *If `G` is a connected
cubic graph and `G ≠ K₄`, then* `Z(G) ≤ α(G) + 1`.

`Z(G)`, the **zero forcing number**, is the least size of a set `S` of vertices
such that colouring `S` and repeatedly applying "a coloured vertex with exactly
one uncoloured neighbour colours that neighbour" eventually colours everything.
`α(G)` is the largest independent set. Theorem 18 of that paper reduces subcubic
graphs to cubic ones, so cubic graphs are the whole conjecture.

Known: proved for claw-free cubic graphs and for orientably one-face-embeddable
cubic graphs; `Z ≤ α + 2` asymptotically almost surely; an infinite family
attains equality. **No exhaustive computational verification is reported.**

## Why this target rather than another

The inventory is below. Two reasons decided it.

**`Z` and `α` are integers.** The two worst recurring failures in this
repository — the `n = 14` star and the `K₃,₃` scan — were both floating-point
equality artifacts, on conjectures that are exactly tight on an infinite family.
This conjecture is *also* exactly tight on an infinite family, and that now
costs nothing.

**The refutation literature structurally avoids it.** Roucairol and Cazenave
state that they address only the conjectures "that do not involve solving an
NP-hard problem". `Z` and `α` are both NP-hard, so the deep-RL and Monte Carlo
work that has swept the spectral conjectures has never been pointed here.

## Instrument validation

`zero_forcing.py` reproduces ten published values — `Z(P₅)=1`, `Z(C₆)=2`,
`Z(K₅)=4`, `Z(K₄)=3`, `Z(K₃,₃)=4`, `Z(K₂,₄)=4`, `Z(Petersen)=5`, `Z(Q₃)=4`,
`Z(prism)=3`, `Z(K₁,₅)=4`, with the matching `α`. The three that matter are the
ones a wrong implementation would get wrong quietly:

* **K₄** must come out as a violation (`Z = 3 > α + 1 = 2`) — it is the one
  graph the conjecture excludes;
* **K₃,₃** and **Petersen** must come out at equality (`Z = α + 1`).

All three behave correctly.

## Generating the graphs: the count check earned its keep immediately

McKay's public collections contain no complete list of connected cubic graphs,
so `cubic_graphs.py` generates them: every connected cubic multigraph on `n+2`
vertices arises from one on `n` by subdividing two edge slots and joining the
new vertices, starting from the theta graph and the dumbbell.

The rule was *assumed*. It is checked against OEIS A002851 at every level, and
the run aborts on the first disagreement — because a silently incomplete
enumeration turns "no counterexample found" into a false statement, and is the
one failure a scan cannot detect from the inside.

**The first version failed that check on its third level**: keeping only simple
graphs produced 4 of the 5 connected cubic graphs on 8 vertices. One of them is
reachable only through a multigraph on 6 vertices. With multigraph intermediates
the counts match exactly: 1, 2, 5, 19, 85 for `n = 4…12`.

## Result so far

**Every connected cubic graph on ≤ 12 vertices** (112 graphs): `K₄` is the only
violation, exactly as the conjecture asserts. Larger orders are running.

The graphs attaining equality, `Z = α + 1`:

| n | α | structure |
|---|---|---|
| 6 | 2 | prism `CL₃` |
| 6 | 3 | `K₃,₃` |
| 8 | 3 | three graphs, one of them the Möbius ladder `M₈` |
| 10 | 4 | three graphs, one of them the **Petersen graph** |
| 12 | 4, 5 | three graphs |

## What went wrong (continued)

**9. A structural characterisation was inferred from a pattern, asserted as if
classical, and built on — and it was false.** The tight cases at `n = 12` all had
`α = 4 = n/3`, which is the minimum for cubic graphs of that order. From this
came the claim that `α = n/3` forces the vertex set to partition into triangles,
so that the frontier of the conjecture is exactly the *triangle-replaced* graphs
`truncate(H)` for `H` cubic on `n/3` vertices. That is an attractive reduction —
it turns `n` into `n/3` and reaches `n = 30` in a second — and `triangle_replaced.py`
was written around it.

It found **zero** tight cases, contradicting the exhaustive scan, which had
found three at `n = 12`. The implementation was correct: `truncate(K₄)` is the
truncated tetrahedron, `α = 4` — but `Z = 4`, not 5, so it is not tight at all.
Both halves of the inference were wrong. Equality does not require minimum `α`
(there is a tight graph at `n = 12` with `α = 5`), and `α = n/3` does not force a
triangle partition.

The error was caught only because exhaustive ground truth existed to contradict
it. Had the frontier family been attacked first — as its cube-root efficiency
made tempting — it would have returned "no counterexample in a family reaching
n = 30" and that would have been reported as progress on a family that does not
contain the tight cases at all.

The lesson is not "check your algebra". It is that a *reduction of the search
space* is a mathematical claim, and carries the same burden of proof as a
result. Speed makes them attractive precisely in proportion to how much they
throw away.

**10. The search optimised against the wrong adversary, then built a wall around
the goal.** `cubic_search.py` explores orders the exhaustive scan cannot reach,
and needed three objectives before one worked.

`Z − α` takes about three distinct values over the whole search space, so it is
flat almost everywhere and an annealer on it random-walks. This is the mirror
image of the Laplacian problem: there the quantities were continuous and
rounding *manufactured* counterexamples; here they are integers, nothing is ever
manufactured, and the difficulty is that nothing moves.

*First objective*: count how many random `(α+1)`-subsets force. It reached zero
within seconds — and the greedy construction then found a forcing set anyway,
every time. Random subsets are a weak adversary, so the search had been
optimising against the wrong opponent and its apparent progress meant nothing.

*Second objective*: add greedy as a penalty, applied only when the random count
hit zero. This built a cliff. States with no random hits were punished by
`+10000`, so the annealer learned to avoid exactly the region worth exploring,
and parked at "one random subset forces" indefinitely. **A hard penalty on the
goal state is a barrier around the goal.**

*Third objective*: count greedy successes as an ordinary graded term,
`1000·(greedy hits) + 10·(random hits) + mean reach`. Every layer now moves, and
the search descends through the greedy dimension instead of avoiding it.

---

# Third target: the TxGraffiti product conjecture — a result

Started 2026-07-26, after dropping the cubic scan.

## The conjecture

**Conjecture** (TxGraffiti; stated open in [arXiv:2507.17780](https://arxiv.org/abs/2507.17780)
and [arXiv:2409.19379](https://arxiv.org/abs/2409.19379)). *For connected `G, H`
with at least two vertices,* `γ_t(G □ H) ≥ γ(G × H)`, *where* `□` *is the
Cartesian product and* `×` *the direct product.*

Chosen because a search for work citing it turns up **nothing** — no proofs, no
partial results, no counterexample hunts. Compare the Laplacian bounds, where
two groups had swept the field and a third closed it six weeks before this
project started.

## Result: the conjecture holds with equality for H = K₂, for every graph G

> **Theorem.** `γ_t(G □ K₂) = γ(G × K₂)` for every graph `G`.

Both products have vertex set `V(G) × {0,1}`. In the prism `G □ K₂`,
`(v,i) ~ (v,1−i)` and `(u,i) ~ (v,i)` when `u ~ v`. In the bipartite double
cover `G × K₂`, `(u,i) ~ (v,1−i)` when `u ~ v`, and nothing else. Let
`φ(T) = {(v,1−i) : (v,i) ∈ T}` — a size-preserving involution on subsets. Then

```
T totally dominates G □ K₂   ⟺   ∀(v,i):  (v,1−i) ∈ T   or  (u,i) ∈ T for some u ~ v
φ(T) dominates G × K₂        ⟺   ∀(v,i):  (v,i) ∈ φ(T)  or  (u,1−i) ∈ φ(T) for some u ~ v
```

and by the definition of `φ`, `(v,i) ∈ φ(T) ⟺ (v,1−i) ∈ T` and
`(u,1−i) ∈ φ(T) ⟺ (u,i) ∈ T`. **The two conditions are the same condition.** So
`φ` is a size-preserving bijection between the total dominating sets of the
prism and the dominating sets of the double cover, and the minima coincide. ∎

`prism_identity.py` checks the bijection *on subsets*, not on the two numbers:
for eight graphs including odd cycles, the Petersen graph and a disconnected
one, every subset up to size 4 gives the same answer to both questions. Zero
mismatches.

### How it was found

Not by searching for it. The exhaustive scan over pairs reported which pairs are
**tight**, and `K₂` appeared in more tight pairs than anything else. Fixing
`H = K₂` and varying `G` gave equality for **100%** of the graphs tried — all
995 connected graphs on ≤ 7 vertices and all 182 structured and random graphs on
8–20 vertices. A conjectured *inequality* that is an *equality* everywhere is
not a coincidence, and the proof followed from writing down the two adjacency
conditions side by side.

### Relation to what is known

Azarija, Henning and Klavžar (*(Total) Domination in Prisms*, Electron. J.
Combin. 24 (2017)) prove `γ_t(G □ K₂) = 2γ(G)` for **bipartite** `G`, using
hypergraph transversals, and note that bipartiteness is essential: for every
`k ≥ 1` some non-bipartite `G` has `γ_t(G □ K₂) = 2γ(G) − k`.

The identity contains both statements. For bipartite `G` the double cover splits
into two copies of `G`, so `γ(G × K₂) = 2γ(G)` and their theorem follows. For
non-bipartite `G` the double cover is connected, and `γ(G × K₂)` *is* the smaller
quantity their remark points at. The smallest case is `C₇`: `γ(C₇) = 3` so
`2γ = 6`, but `C₇ × K₂ = C₁₄` and `γ(C₁₄) = 5`. Exactly 11 connected graphs on
≤ 7 vertices have `γ(G × K₂) < 2γ(G)`, all non-bipartite, all with deficit 1.

The proof needs no hypergraph machinery and **no hypothesis on `G` at all** —
not connectedness, not bipartiteness.

### The literature check, done properly this time

The two papers that own this corner of the subject are Azarija–Henning–Klavžar
(2017) and the Goddard–Henning note answering it (*A Note on Domination and
Total Domination in Prisms*). Both were read, not just searched for.

* **Azarija–Henning–Klavžar** prove the bipartite case by hypergraph
  transversals. The paper does not mention the direct product or double covers
  anywhere.
* **Goddard–Henning** give a short proof of the same bipartite case, and close
  it with this remark:

  > "One can also establish the above theorem by noting that when `G` is
  > bipartite, the open neighborhood hypergraph of the prism `G □ K₂` is
  > isomorphic to two disjoint copies of the closed neighborhood hypergraph of
  > `G`."

That remark is the identity above, **stated only in the bipartite case**. The
general form is: the open neighbourhood hypergraph of `G □ K₂` is isomorphic to
the closed neighbourhood hypergraph of `G × K₂`, for every `G`. When `G` is
bipartite the double cover splits into two copies of `G` and their sentence
falls out.

So the honest position is:

* the identity is **not stated** in either of the two papers that would state
  it, nor anywhere the searches reached;
* it is **one step** from a published closing remark, and a specialist reading
  that sentence would produce it in a minute;
* it is therefore best treated as folklore-adjacent — worth recording for what
  it settles, not worth claiming as a discovery.

What it does settle is the TxGraffiti conjecture for `H = K₂`, which is a
statement nobody appears to have made.

### Why the same trick does not generalise

To get `γ_t(G □ H) ≥ γ(G × H)` from a map `ψ` on vertices, it is enough that for
every `x ∈ V(G × H)` there is a `y ∈ V(G □ H)` with
`ψ(N_{G□H}(y)) ⊆ N_{G×H}[x]`. Taking `ψ = id` and `x = (v,h)`, the candidate
`y = (v,h'')` with `h'' ~ h` works for the `{(u,h'') : u ~ v}` part but forces
`N_H(h'') = {h}` for the rest — i.e. `h` must have a pendant neighbour in `H`.
`H = K₂` is exactly the case where that holds for every `h`.

Pushing on this gives a clean statement of why the approach stops there. The
containment argument works at `(v,h)` as soon as **`v` has a pendant neighbour in
`G` or `h` has one in `H`** — take `y = (v', h)` for a pendant neighbour `v'` of
`v`, or `y = (v, h')` for a pendant neighbour `h'` of `h`, and in either case
`N_{G□H}(y) ⊆ N_{G×H}[(v,h)]`. Whenever this holds at every pair, *every* total
dominating set of `G □ H` is already a dominating set of `G × H`, with no
modification at all. Verified: over four choices of `H` and seven of `G`, every
total dominating set of size ≤ 6 also dominates the direct product. Zero
violations.

But the hypothesis is nearly vacuous. If some `v` has no pendant neighbour then
*every* `h` must have one; a degree-1 vertex `u` has a single neighbour `w`, and
for `u` to have a pendant neighbour `w` must have degree 1 too, so `{u,w}` is a
`K₂` component. **For connected `H` the condition holds only for `H = K₂`** — and
symmetrically for `G`. So this is not a partial result that might be widened; it
is a proof that the containment approach reaches exactly `K₂` and stops.

The general conjecture needs a different idea, not a better `ψ`.

## Exhaustive status of the general conjecture

| range | pairs | result |
|---|---|---|
| both graphs on ≤ 5 vertices | 465 | no counterexample; 29% exactly tight |
| both graphs on ≤ 6 vertices | 10,153 | no counterexample; 7.5% exactly tight |
| **all pairs with `\|G\|·\|H\| ≤ 42`** (i.e. everything but `7×7`) | **131,279** | **no counterexample** |
| the `7×7` pairs | 364,231 | running |
| local search, 22 shapes up to products of 48 vertices | — | no counterexample |
| `H` a path or cycle of order 3–8, `G` to 11 vertices | 1,144 | no counterexample |

Minimum slack is `0` throughout: the conjecture is tight but never violated.

### Where the conjecture is tight

Tightness falls off sharply with size, which is why a counterexample — if there
is one — should be expected small:

| `(|G|,|H|)` | tight | | `(|G|,|H|)` | tight |
|---|---|---|---|---|
| (2, *) | **100%** | | (4,6) | 16.4% |
| (3,4) | 100% | | (5,5) | 6.5% |
| (3,5) | 76.2% | | (5,6) | 4.5% |
| (3,6) | 54.5% | | (6,6) | 2.9% |
| (4,4) | 42.9% | | | |

The `(2, *)` row is the `H = K₂` theorem.

Controlling for size (both graphs on ≥ 5 vertices, so the small-graph effect is
removed), the properties that predict tightness are:

| property | share of tight pairs | share of all pairs | lift |
|---|---|---|---|
| both are trees | 3.0% | 0.5% | **5.9×** |
| either is a tree | 45.8% | 13.0% | **3.5×** |
| both bipartite | 8.3% | 2.8% | 2.9× |
| either bipartite | 72.4% | 30.2% | 2.4× |
| either has a pendant vertex | 84.4% | 70.5% | 1.2× |

And conditioning on one side lying in a named class, with both sides on ≥ 4
vertices:

| class of one side | tight rate |
|---|---|
| **cycle** | **38.6%** |
| **path** | **33.8%** |
| other tree | 8.0% |
| star | 6.4% |
| complete | 4.5% |
| other | 2.9% |

**This contradicts the intuition the proof of the `K₂` case suggests.** The
containment argument works wherever a vertex has a pendant neighbour, which
points at stars and pendant-rich graphs. The data points somewhere else
entirely: stars are barely above baseline (6.4%), and the tight pairs
concentrate on paths and cycles. `K_{1,k}` is tight for 100% of partners at
`k = 1` (that is `K₂`), 56% at `k = 2` (that is `P₃`), and then **6.4% flat**
for every `k ≥ 3` — so the star family is not the continuation of the `K₂` case
at all. Whatever makes `K₂` work does not extend along stars.

The search is now prioritised on paths and cycles. That is a prioritisation, not
a claim that counterexamples must live there — see failure 9 for what happens
when the two are confused.

### Tightness along paths and cycles does not decay

Fixing `H` and running `G` over all connected graphs on 4–6 vertices plus
structured and random graphs up to 11:

| `H` | tight | | `H` | tight |
|---|---|---|---|---|
| `P₃` | 60.1% | | `C₃` | 60.1% |
| `P₄` | 39.5% | | `C₄` | 51.9% |
| `P₅` | 26.8% | | `C₅` | 30.9% |
| `P₆` | 38.1% | | `C₆` | 37.4% |
| `P₇` | 37.0% | | `C₇` | 29.6% |
| `P₈` | 33.3% | | `C₈` | 29.6% |

No counterexample; minimum slack `0` throughout.

The rate **flattens out around 30–40% and stays there** as `|H|` grows. Contrast
the stars, where it collapses to 6.4% at `k = 3` and never recovers. So paths and
cycles are not merely where tightness happens to be common at small sizes — the
conjecture stays on its boundary along these two families indefinitely, which is
what makes them the place to look and the place a proof has to survive.

### A source that turned out not to be one: the Optimist

Davila's *Optimist* ([arXiv:2411.09158](https://arxiv.org/abs/2411.09158),
[repository](https://github.com/RandyRDavila/The-Optimist)) looked promising as a
successor vein — a newer automated conjecturer, less picked over than Graffiti.
It is not a source of open problems.

The repository is a single demo notebook. Its output is **raw, unfiltered
generated conjectures**, not a curated list. Extracting all 103 distinct ones and
running them over the 995 connected graphs on ≤ 7 vertices:

* **80 of 103 are refuted immediately** by that scan;
* the 23 survivors are, on inspection, mostly known theorems — `α ≤ n − μ`
  (Gallai), `α = n − μ` for bipartite graphs (König), `α ≥ γ` (a maximal
  independent set dominates), `α ≥ n/2` for bipartite graphs.

Which is exactly what the paper says the system does: "rediscover established
theorems and propose novel inequalities". Refuting the rest is no contribution —
the system would refute them itself given a larger database.

The lesson for target selection: **what made Written on the Wall II unusually
good was that it was curated, formalised, and obscure — all three.** Raw
generator output fails the first, AutoGraphiX fails the second, and the famous
conjectures fail the third.

---

# Fourth target: the Graffiti.pc conjectures, from their formalisation

## A better version of a source that was rejected earlier

The Graffiti corpus was ruled out early in this project: Roucairol and Cazenave
document that its definitions disagree between sources, and that they produced
two false refutations from mis-defined invariants. Definition risk, not
difficulty, was the reason to stay away.

Google DeepMind's [`formal-conjectures`](https://github.com/google-deepmind/formal-conjectures)
removes that risk. It holds **631 open conjectures formalised in Lean**, of which
178 are combinatorics, and its `WrittenOnTheWallII` directory contains **22 of
DeLaVina's Graffiti.pc conjectures tagged `category research open`**. The
definition *is* the Lean definition; there is nothing left to misread.

By source: Erdős problems 353, Wikipedia 123, Green's problems 47, papers 24,
**Written on the Wall II 22**, OEIS 17, MathOverflow 11, arXiv 10, rest 24.

## Choosing among the 22 by how worked-over each theme already is

Counting mentions on DeLaVina's own page of *resolved* conjectures:

| theme | mentions | conjectures |
|---|---|---|
| induced bipartite `b(G)` | 40 | 19, 40, 59, 198 |
| Hamiltonicity | 25 | 194, 198, 200, 217 |
| forest number `f(G)` | 24 | 40, 59, 61, 65 |
| local independence `l(v)` | 22 | 2, 100, 141, 160, 194 |
| total domination | 4 | 291 |
| **induced tree** | **2** | 141, 142, 144, 145, 146 |
| **residue** | **0** | 59, 61, 217, 291 |
| **well totally dominated** | **0** | 314, 316, 322 |

The bottom three themes are where nothing has been resolved, so that is where to
look. This is target selection done before writing code rather than after — the
step whose absence cost the whole Laplacian result.

## Result: conjecture 322 is true, and its hypothesis is nearly vacuous

> **WOWII 322.** *Let `G` be connected on `n ≥ 5` vertices. If `max_v l(v) ≤ 1`,
> where `l(v) = α(G[N(v)])`, then `G` is well totally dominated.*

**Theorem.** True, and `n ≥ 5` is unnecessary — it holds for every connected
graph on at least two vertices.

*Step 1: the hypothesis forces `G` to be complete.* `l(v) ≤ 1` says `N(v)` is a
clique. If `u ~ v ~ w` with `u ≠ w`, then `u, w ∈ N(v)`, so `u ~ w`. Adjacency is
transitive, every component is a clique, and connectedness gives `G = K_n`.

*Step 2: `K_n` is well totally dominated.* Every pair `{u,v}` totally dominates;
a singleton does not, since `u` has no neighbour in `{u}`; and any set of size
`≥ 3` properly contains a dominating pair, so is not minimal. All minimal total
dominating sets have size exactly 2. ∎

Worth almost nothing as mathematics — two lines, unwritten presumably because
nobody looked. What it illustrates is a property of machine-generated
conjectures: `max_v l(v) ≤ 1` *reads* as a mild local condition and is in fact
the strongest possible one, collapsing the statement to a single family. Of the
995 connected graphs on ≤ 7 vertices, exactly **6** satisfy it: `K₂` through
`K₇`.

Checking what a hypothesis actually admits, before spending a search budget on
it, is cheap.

## Result: conjecture 141 is true — in its original, stronger form

> **Correction.** What was proved first was the statement *as formalised in
> Lean*, `tree(G) ≥ ⌊girth/2⌋ − 1 + max_v l(v)`. DeLaViña's original asks for
> `(1/2)·girth`, which for **odd** girth is stronger by 1 after rounding — on
> `C₇` the old argument gave 4 where the original needs 4.5, i.e. 5. Two gaps
> had to be closed, at girth 5 and at odd girth ≥ 7; both are, and the result
> below now covers the original. The gap was invisible to computation: exactly
> **two** of the 12,112 connected graphs on ≤ 8 vertices have odd girth ≥ 7, so
> no exhaustive check in range can tell the two readings apart. See
> `wowii141.py` for the corrected proof — the extra ingredient is counting the
> ball more carefully, plus one vertex at distance `r+1` in the tight case.

### The original proof (of the weaker, formalised reading)

> **WOWII 141.** *`tree(G) ≥ ⌊girth(G)/2⌋ − 1 + max_v l(v)`, with `tree(G)` the
> order of a largest induced tree and `l(v) = α(G[N(v)])`.*

**Theorem.** True. Two constructions, one per range of the girth, and the ranges
meet exactly.

> **Lemma A.** `tree(G) ≥ 1 + max_v l(v)`.
>
> *Proof.* Pick `v` attaining the maximum and let `S` be a maximum independent
> set of `G[N(v)]`. In `G[{v} ∪ S]` the only edges join `v` to `S`, since `S` is
> independent — a **star centred at `v`**, a tree of order `1 + l(v)`. ∎

*Case `girth ≤ 5`* (including acyclic, girth `0` by the convention): then
`⌊g/2⌋ − 1 ≤ 1`, and Lemma A already gives it.

*Case `girth g ≥ 6`*: no triangles, so `N(v)` is independent for every `v`,
hence `l(v) = deg(v)` and `max_v l(v) = Δ`. Put `r = ⌊g/2⌋ − 1 ≥ 2` and take `v`
of degree `Δ`.

* **`B(v,r)` induces a tree.** A cycle inside a ball of radius `r` has length at
  most `2r+1`: it contains an edge `xy` outside the BFS tree from `v`, whose
  fundamental cycle has length `≤ dist(v,x) + dist(v,y) + 1 ≤ 2r+1`. Here
  `2r+1 = 2⌊g/2⌋ − 1`, which is `g−1` for even `g` and `g−2` for odd `g` — below
  the girth either way, so there is no cycle.
* **`|B(v,r)| ≥ Δ + r`.** The ball induces a tree while `G` has a cycle, so
  `B(v,r) ≠ V`, hence `ecc(v) > r` and some vertex sits at each distance
  `1, …, r`. Counting `v`, its `Δ` neighbours and one vertex at each distance
  `2, …, r` gives `1 + Δ + (r−1) = Δ + r`.

So `tree(G) ≥ Δ + ⌊g/2⌋ − 1`. ∎

Both ingredients are textbook — the neighbourhood star, and balls of small radius
being locally tree-like. What had not been done is pointing them at this
statement. **The split point is forced**: Lemma A works exactly when
`⌊g/2⌋ − 1 ≤ 1`, i.e. `g ≤ 5`, and the ball needs `r ≥ 2`, i.e. `g ≥ 6`. Nothing
lies between. Verified on all 995 graphs on `≤ 7` vertices and on Heawood,
Möbius–Kantor, Pappus, Desargues, `C₇`, `C₉`, `C₁₂`, `C₂₀`, grids and `Q₄`.

## Result: conjecture 65 is true, trivially

> **WOWII 65.** *`distMin(A) + ⌈distMin(M)/3⌉ ≤ f(G)`, with `A` the
> minimum-degree vertices and `M` the maximum-degree ones.*

`distMin(S)` is the minimum over `v ∉ S` of `dist(v,S)`. In a **connected** graph
any nonempty proper `S` has a vertex adjacent to it, so `distMin(S) = 1`; and
`distMin(V) = 0`. It takes no other value. So the left-hand side is at most
`1 + ⌈1/3⌉ = 2`, while `f(G) ≥ 2` for any graph with an edge. ∎ Over the 995
graphs the left-hand side takes only the values `{0, 2}`.

## Result: conjecture 316 is true

> **WOWII 316.** *If `avgdeg(Ḡ) ≤ |P|` with `P` the pendant vertices, then `G` is
> well totally dominated.*

*Step 1.* With `p = |P|`, core size `r = n − p` and `m` edges,
`avgdeg(Ḡ) = (n−1) − 2m/n` makes the hypothesis `m ≥ (r+p)(r−1)/2`. Each pendant
contributes one edge and the core spans at most `C(r,2)`, so
`m ≤ p + r(r−1)/2`. Combining gives `p ≥ (r−1)p/2`, hence **`r ≤ 3`** whenever
`p > 0`, and each case pins the core: `r = 1` is the star; `r = 2` forces the
core edge, since otherwise each core vertex takes its pendants into a separate
component; `r = 3` needs `m ≥ 3 + p` and so all three core edges. With `p = 0`
the hypothesis is `m ≥ n(n−1)/2`, i.e. `G = K_n`. So it admits exactly the
complete graphs and `K₁`, `K₂`, `K₃` with pendants.

*Step 2.* Each of those is well totally dominated — `K_n` as in 322, and for a
`K₃` core all three carrying pendants gives the unique minimal set `{u,v,w}` of
size 3, while otherwise a forced vertex plus a neighbour is minimal of size 2. ∎

## Result: conjecture 61 settled for diameter ≤ 3

> **WOWII 61.** *`f(G) ≥ residue(G) + ⌈diam(G)/3⌉`, with `f` the order of a
> largest induced forest.*

> **Lemma.** `f(G) ≥ α(G) + 1` for every connected graph on `≥ 2` vertices.
>
> *Proof.* Let `I` be a maximum independent set and `v` a vertex outside it. In
> `G[I ∪ {v}]` the only edges run between `v` and `N(v) ∩ I` — a star centred at
> `v` plus isolated vertices, which is a forest. ∎

**Corollary.** Conjecture 61 holds whenever `diam(G) ≤ 3`, since then
`⌈diam/3⌉ = 1` and Favaron–Mahéo–Saclé give `residue ≤ α`.

With the second free bound `f ≥ diam + 1` these settle everything in reach: 893
of the 995 graphs have diameter `≤ 3`, all but 5 of the remaining 102 fall to
`f ≥ diam + 1`, and those 5 hold anyway. **Diameter `≥ 4` is open.** The route
would be a strengthening: a set `S` of vertices pairwise at distance `≥ 3` and
disjoint from a maximum independent set `I` has no two members sharing a
neighbour, so `G[I ∪ S]` is a disjoint union of stars and `f ≥ α + |S|`; a
shortest path of length `d` supplies `⌊d/3⌋ + 1` such vertices, and making them
disjoint from *some* maximum independent set is the gap.

The α-strengthened form `f ≥ α + ⌈diam/3⌉` is **false** — 8 violations, smallest
on 6 vertices — so `residue` is doing real work and no proof can route through
`α`.

## The free-bound toolkit, and what it leaves

Every result above came from the same move: find a lower bound on the invariant
that costs nothing, then see which range of the conjecture it already covers.
All verified over the 995 graphs with zero violations:

| bound | reason |
|---|---|
| `f(G) ≥ α + 1` | a maximum independent set plus one vertex induces a **star** |
| `f(G) ≥ diam + 1` | a shortest path is **induced** |
| `tree(G) ≥ diam + 1` | a shortest path is an induced **tree** |
| `tree(G) ≥ 1 + max_v l(v)` | `v` plus a maximum independent set of `N(v)` is a star |
| `tree(G) ≥ girth − 1` | a shortest cycle is **chordless**; drop one vertex |
| `b(G) ≥ f(G)` | a forest is bipartite |
| `path(G) ≥ diam + 1` | a shortest path is an induced path |

Applied to the three remaining induced-tree conjectures:

| conj | covered by a free bound | residual |
|---|---|---|
| 142 `(2/3)girth + eccSet(B) ≤ tree` | 985/995 | **10** |
| 144 `girth − 1 + ecc(centre) ≤ tree` | 987/995 | **8** |
| 146 `2·eccSet(B) ≤ tree·rad(G²)` | 982/995 | **13** |

Every residual graph satisfies its conjecture; none is refuted. They sit at
girth 4–6 and are mostly tight, which is where the free bounds run out.

### A new free bound: hanging a path off a shortest cycle

Pushing on `#144` produced one more bound, and it needs a hypothesis that the
first attempt missed.

> **Lemma.** If `girth(G) = g ≥ 5` and `C` is a shortest cycle, then
> `tree(G) ≥ (g − 1) + ecc(C)`, where `ecc(C) = max_w dist(w, C)`.
>
> *Proof.* A shortest cycle is chordless. Take `w` with `dist(w,C) = ecc(C)` and
> a shortest path `Q = w = q_d, …, q_1, q_0 = c ∈ C`. For `i ≥ 2`,
> `dist(q_i, C) ≥ 2`, so `q_i` has no neighbour on `C`. And **`q_1` has at most
> one neighbour on `C`**: two neighbours `c_i, c_j` would close a cycle through
> `q_1` and the shorter arc, of length `≤ ⌊g/2⌋ + 2`, which is `< g` exactly when
> `⌈g/2⌉ > 2`, i.e. `g ≥ 5`. So `C ∪ Q` induces a cycle with one pendant path;
> deleting any vertex of `C` other than `c` leaves an induced tree of order
> `(g − 1) + d`. ∎

**The hypothesis is `girth ≥ 4`, not `≥ 5`.** The first version of the argument
stopped at 5, which was too cautious. Let `N = N(q_1) ∩ C`. For `g ≥ 4`, `N` is
independent in `C` — two adjacent members would give a triangle — and the arc
argument gives `|N| ≤ 1` for `g ≥ 5` and `|N| ≤ 2` for `g = 4`, the two being
diagonal. Both cases build the tree, with different deletions:

* `|N| = 1`: delete a vertex of `C` **outside** `N`. The remaining `g−1` vertices
  form a path and `q_1` attaches at the single vertex of `N`.
* `|N| = 2` (so `g = 4`): delete **one member of `N`**. The remaining three
  vertices form a path and `q_1` attaches only at the other member.

At girth 3 it genuinely fails, and `K₄` shows why: there `q_1` is adjacent to all
three cycle vertices, since a triangle no longer contradicts the girth. Verified
over the 995 graphs — **28 violations at girth 3, zero at girth 4 or above.**

This takes `#144`'s residual from 8 graphs to **3**, using only proved bounds.

**What blocks the last three is exactly `ecc(C)` versus `ecc(centre)`.** The
lemma delivers distance to the shortest cycle; the conjecture asks for distance
to the centre, and these differ on 7 of the 971 graphs with a cycle. In the
smallest residual case — `C₅` with a pendant on each of two adjacent cycle
vertices — `ecc(C) = 1` while `ecc(centre) = 2`, and the induced tree of order 6
that does exist is the cycle minus a vertex **plus both pendants**. So the
missing ingredient is hanging paths at *several* cycle vertices at once, not one;
the single-path lemma cannot see it.

### The caterpillar bound, and conjecture 144 closed in range

The single-path lemma misses cases where several paths hang off the cycle at
once. Generalising costs nothing:

> **Lemma (caterpillar).** Let `C` be a shortest cycle with `girth ≥ 4`, let
> `x ∈ C`, and let `P = C − x` (an induced path on `g−1` vertices). If `A` is a
> set of vertices outside `C` that is **independent** and in which every member
> has **exactly one** neighbour in `P`, then `G[P ∪ A]` is an induced tree — a
> caterpillar — of order `(g−1) + |A|`. ∎

Verified with **zero violations** over the 995 graphs. Adding it:

| conj | residual before | after |
|---|---|---|
| 142 | 6 | **2** |
| 144 | 3 | **0** |
| 146 | 13 | 11 |

**Conjecture 144 is settled over every connected graph on ≤ 7 vertices by proved
lower bounds alone**, and over all 12,111 connected graphs on ≤ 8 vertices the
three cheap bounds leave 29, of which the caterpillar closes all but **4**.

The smallest case it was built for makes the point: `C₅` with a pendant on each
of two adjacent cycle vertices has `ecc(C) = 1` but `ecc(centre) = 2`, and the
induced tree of order 6 is the cycle minus a vertex **plus both pendants** — two
attachments, which is exactly what the single-path version cannot see.

This is not yet a proof of `#144`. What is proved is the bound; what is checked
is that it beats the conjecture's right-hand side on everything up to 8
vertices, bar four graphs. Closing it needs `|A| ≥ ecc(centre)` for some choice
of `C` and `x`, which is where the argument now sits.

### Conjecture 133 holds for every graph containing a 4-cycle

> **WOWII 133.** *`path(G) ≥ rad(G) + ⌊l_avg(G)⌋^{c}` where `c = 0` if `G`
> contains a 4-cycle and `1` otherwise.*

The exponent does all the work. **If `G` has a 4-cycle then `c = 0`, the power is
`1`, and the claim is `path(G) ≥ rad + 1`** — which the free bound
`path ≥ diam + 1 ≥ rad + 1` gives immediately. That is 905 of the 995 graphs,
settled outright.

The remaining 90 graphs have no 4-cycle and face `path ≥ rad + ⌊l_avg⌋`. The same
free bound settles 86 of them; the residual is **one** graph, `n = 7` with
`rad = diam = 2`, `⌊l_avg⌋ = 2` and `path = 5` against a right-hand side of `4`.

### Conjecture 40 holds for every graph with a Hamiltonian path

> **WOWII 40.** *`f(G) ≥ ⌈(p(G) + b(G) + 1)/2⌉`, with `p` the path cover number
> and `b` the largest induced bipartite subgraph.*

Since `b ≥ f` always (a forest is bipartite), the conjecture is equivalent to
the cleaner

```
2 f(G)  ≥  p(G) + b(G) + 1.
```

**When `p(G) = 1` — that is, when `G` has a Hamiltonian path — the free bounds
`f ≥ α + 1` and `f ≥ diam + 1` settle it, all 851 of 851 such graphs.** The
residual is **3 graphs**, all with `p = 2`, `b = f = 6`, where the conjecture
holds comfortably (`12 ≥ 9`) but the free bounds fall one short.

Worth noting in passing: the conjecture implies `f ≥ p + 1`, which is itself
true with no violations — and `b − f ∈ {0,1,2}` across all 995 graphs, so the
gap between the two invariants is never large.

### Two more free bounds, and conjecture 19 closed in range

Both come from the same observation as Lemma A, pushed a little.

> **Lemma B.** `b(G) ≥ α(G) + α(G − I)` for a maximum independent set `I`.
>
> *Proof.* `I` and a maximum independent set `J` of `G − I` are disjoint and each
> independent, so `G[I ∪ J]` is bipartite. ∎

> **Lemma C.** `f(G) ≥ α(G) + |A|` whenever `A` is an independent set of vertices
> outside a maximum independent set `I`, each with **exactly one** neighbour in
> `I`.
>
> *Proof.* `G[I ∪ A]` has edges only from `A` to `I`, one per member of `A`, so it
> is a disjoint union of stars — a forest. ∎ (Lemma A is the case `|A| = 1`.)

Both verified with zero violations over the 995 graphs. Their effect:

| conj | residual before | after |
|---|---|---|
| **19** `⌊avg ecc + max l(v)⌋ ≤ b` | 148 | **0** |
| **59** `f ≥ ⌈√(residue·b)⌉` | 74 | **22** |

**Conjecture 19 is settled over every connected graph on ≤ 7 vertices by proved
lower bounds alone** — and it had the largest residual of the whole family
before Lemma B.

### More free bounds, and where each conjecture stands

Four more, each verified with zero violations over the 995 graphs:

| bound | reason |
|---|---|
| `Ls(G) ≥ Δ` | a spanning tree keeping every edge at a maximum-degree vertex |
| `b(G) ≥ α + 1` | the star of the previous lemma is bipartite |
| `path(G) ≥ girth − 1` | a chordless shortest cycle minus a vertex is a path |
| `γ_t(G) ≤ n − Δ + 1` | a maximum-degree vertex and one neighbour cover its closed neighbourhood |

Coverage of the remaining conjectures by free bounds alone:

| conj | covered | residual |
|---|---|---|
| 133 | 994/995 | **1** |
| 40 | 992/995 | **3** |
| 2 | 986/995 | 9 |
| 291 | 989/994 | 5 |
| 59 | 921/995 | 74 |
| 19 | 847/995 | 148 |

`#133` is one graph away and `#40` three. In all nine residual cases of `#2` the
graph has `Ls = Δ + 1` exactly — but `Ls ≥ Δ + 1` is false in general (stars and
cycles give `Ls = Δ`), so that observation does not lift to a bound.

One observation from the failed attempt on 144: **`ecc(centre) = diam − radius`
on 989 of the 995 graphs**, and `diam − radius + 1` on the other 6. That
near-identity suggested `girth ≤ radius + 2` as a sufficient condition; it does
hold on 985 graphs, but the mechanism is not the conjectured one — `diam + 1 ≥ RHS`
fails there. The correlation is real and the proof is not.

## Fifteen of the 22, exhaustively

All connected graphs on ≤ 6 vertices (142 graphs), no counterexample anywhere:

| conj | statement | applicable | holds |
|---|---|---|---|
| 2 | `Ls(G) ≥ 2(l_avg − 1)` | 142 | 142 |
| 40 | `f(G) ≥ ⌈(p(G)+b(G)+1)/2⌉` | 142 | 142 |
| 59 | `f(G) ≥ ⌈√(residue·b)⌉` | 142 | 142 |
| 61 | `f(G) ≥ residue(G) + ⌈diam/3⌉` | 142 | 142 |
| 141 | `tree(G) ≥ ⌊girth/2⌋ − 1 + max_v l(v)` | 142 | 142 |
| 144 | `tree(G) ≥ girth − 1 + ecc(centre)` | 142 | 142 |
| 145 | `2·eccSet(B) ≤ tree(G)·lMin(Ḡ)` | 90 | 90 |
| 194 | `α ≤ 1 + l_avg` ⟹ Hamiltonian path | 84 | 84 |
| 198a | `b(G) ≤ 2 + avg ecc` ⟹ Hamiltonian path | 53 | 53 |
| 200 | `tree(G) = ⌈1+l_avg⌉` ⟹ Hamiltonian path | 38 | 38 |
| 217 | `Ls(G) ≤ 4·[residue=2]+2` ⟹ Hamiltonian path | 85 | 85 |
| 291 | `γ_t ≤ HH-zero-step + freq(min triangles)` | 141 | 141 |
| 314 | triangle-free, induced path ≤ 4 ⟹ WTD | 21 | 21 |
| 316 | `avgdeg(Ḡ) ≤ #pendants` ⟹ WTD | 19 | 19 |
| 322 | `n≥5`, all `l(v)≤1` ⟹ WTD | 2 | 2 |

### Lifting the 7-vertex ceiling

Every exhaustive scan in this repository had been capped at 7 vertices, because
that is where networkx's atlas stops. `allgraphs.py` removes the cap: a graph on
`n` vertices is a graph on `n−1` plus a vertex joined to some subset, so applying
that to a complete list gives a complete list. As with the cubic generator the
rule is checked rather than trusted, against **OEIS A000088** (all graphs) and
**A001349** (connected), with an abort on the first disagreement.

Counts match exactly through `n = 8`: 12,346 graphs of which **11,117 are
connected**, generated in 44 seconds. That is 11× the previous ceiling, and it
applies to every scan here, not just this one.

Reaching it needed one algorithmic change. `has_hamiltonian_path` was a
permutation search — 40,320 per graph at `n = 8`, times 11,117 graphs, which is
hopeless. Replaced with Held–Karp over subsets (`O(2ⁿn²)`), validated against
`P₅`, `C₆`, `K₆`, `K_{1,3}`, `K_{2,4}`, `K_{3,3}` and the Petersen graph.

### Where these conjectures are tight, and why that produced nothing here

The method that produced the `H = K₂` theorem was: read the *equality* cases,
not the verdict. Applied to the seven inequality conjectures, all are tight
somewhere and none is ever violated, with tight rates from 0.5% (`#2`) to 32%
(`#59`).

Three classes come out 100% tight — stars for `#59` and `#61`, complete graphs
for `#40`, `#59` and `#61`. But these are the standard extremal families a
conjecture generator fits its bounds to, and each verifies in two lines: on
`K_n`, `f = 2`, `residue = 1`, `diam = 1`, so `#61`'s right-hand side is
`1 + ⌈1/3⌉ = 2`. Equality on the family the bound was built from is not a
structural discovery.

That is the difference from the product conjecture, where `H = K₂` gave equality
for **every** `G` — an identity over an infinite, unremarkable family. Nothing of
that shape appears here.

### All 22 implemented; twelve of them exhaustive to 8 vertices

With the ceiling lifted, twelve conjectures were rerun over **every connected
graph on ≤ 8 vertices (12,112 graphs)**. All hold:

| conj | applicable | holds | | conj | applicable | holds |
|---|---|---|---|---|---|---|
| 59 | 12,112 | 12,112 | | 200 | 671 | 671 |
| 61 | 12,112 | 12,112 | | 291 | 12,111 | 12,111 |
| 141 | 12,112 | 12,112 | | 314 | 81 | 81 |
| 144 | 12,112 | 12,112 | | 316 | 37 | 37 |
| 145 | 10,860 | 10,860 | | 322 | 4 | 4 |
| 194 | 5,621 | 5,621 | | 198a | 516 | 516 |

The remaining seven (19, 65, 100, 133, 142, 146, 160) were then implemented,
completing all 22. Six hold over the 995 graphs on ≤ 7 vertices. The seventh is
below.

### Conjecture 160, as formalised, is false — with a five-vertex counterexample

The small-witness guard from failure 12 fired on `#160`: 489 violations, the
smallest on 5 vertices. It is not a bug this time.

> **Counterexample.** `G = K₅` minus the edges `02` and `13` — that is, `K₅` with
> a perfect matching removed from four of its vertices.

Every quantity is hand-checkable. Degrees are `3,3,3,3,4`, so vertex 4 sees
everything. `N(4) = {0,1,2,3}` induces exactly the 4-cycle `0–1–2–3–0`, so
`l(4) = α(C₄) = 2`, and every other neighbourhood induces a path on three
vertices, also with independence number 2 — so `max_v l(v) = 2`. The triangles
are `{4,0,1}, {4,1,2}, {4,2,3}, {4,3,0}`, all through vertex 4, giving
`max_v T(v) = 4`. The only 4-set inducing a `C₄` is `{0,1,2,3}`, since every
other one contains vertex 4 and so induces at least 5 edges — `c_{C₄}(G) = 1`.
A spanning tree on 5 vertices has at most 4 leaves and the star at vertex 4
attains it, so `Ls(G) = 4`. Then

```
max_v l(v) + max_v T(v) · c_{C₄}(G)  =  2 + 4·1  =  6  >  4  =  Ls(G).
```

**This refutes the statement as formalised; it does not settle DeLaVina's
conjecture 160,** because the formalisation may not render it faithfully — and
there is direct evidence of that risk two files away. Conjecture 133 defines
`cC4` as an *indicator* (`if hasC4 then 0 else 1`) while 160 uses the same name
for a *count*, and conjecture 100's prose says `diam(Ḡ)` where its statement uses
`degreeL2Norm Ḡ`.

The indicator reading does not rescue it either. All three readings fail:

| reading of `c_{C₄}` | violations (n ≤ 6) | smallest |
|---|---|---|
| count, as formalised | 34 | 5 vertices |
| indicator, `0` if `C₄` (conjecture 133's convention) | 53 | 4 vertices |
| indicator, `1` if `C₄` | 29 | 5 vertices |

Whatever conjecture 160 says, the file does not say it. The useful output is a
bug report against the formalisation, not a mathematical claim — and the right
place for it is an issue on `google-deepmind/formal-conjectures`.

### A full prose-vs-statement audit of the 22

Since two problems had already turned up by reading definitions rather than by
searching, all 22 were audited the same way: extract each docstring and each
Lean statement, and compare which mathematical concepts appear in each.

Four came out flagged; three were artefacts of the matching — the docstrings
write `\mathrm{tree}(G)` and `\operatorname{rad}(G)`, which the pattern missed.
**One is real: conjecture 100.**

### Conjecture 100's documentation and its statement disagree

The doc comment reads `α(G) ≤ ⌈(max_v l(v) + ½·diam(Ḡ))/2⌉`. The Lean statement
uses **`degreeL2Norm Gᶜ`** — the square root of the sum of squared degrees — not
the diameter. These are very different quantities. The statement is what is
formalised, so the statement is what is tested here; but anyone reading only the
prose would implement a different conjecture.

### The four Hamiltonian-path conjectures mostly follow from Chvátal–Erdős

`#194`, `#198a`, `#200` and `#217` all conclude "then `G` has a Hamiltonian
path". Chvátal–Erdős gives one for free whenever `κ(G) ≥ α(G) − 1`. Over the
connected graphs on 3–7 vertices:

| conj | satisfy the hypothesis | covered by Chvátal–Erdős | left over |
|---|---|---|---|
| 194 | 500 | 426 (85%) | 74 |
| 198a | 147 | 117 (80%) | 30 |
| **200** | 135 | **133 (99%)** | **2** |
| 217 | 452 | 396 (88%) | 56 |

and the residue is almost entirely `(α, κ) = (3, 1)` — graphs with a cut vertex
and independence number 3. So each of these four reduces to a small, sharply
described class.

**`#217` splits further, and half of it is trivial.** `Ls(G) ≥ 2` always, so when
`residue(G) ≠ 2` the hypothesis is exactly `Ls(G) = 2`. A connected graph with a
vertex of degree `≥ 3` has a spanning tree in which that vertex keeps all three
edges, hence at least three leaves — so `Ls = 2` forces maximum degree `≤ 2`,
i.e. a path or a cycle, both of which have Hamiltonian paths. The whole content
of `#217` is the case `residue(G) = 2` **and** `Ls(G) ≤ 6`.

### The hypotheses rule out the obvious obstruction

A vertex whose deletion leaves `≥ 3` components makes a Hamiltonian path
impossible outright — a path enters and leaves such a vertex once, so it can
join at most two of the pieces. Any graph like that satisfying one of the four
hypotheses would be an immediate counterexample.

Among the 994 connected graphs on 3–7 vertices, **93 have such a vertex**, and
**none of them satisfies any of the four hypotheses**. The Chvátal–Erdős
residual is correspondingly clean: every graph in it has maximum
deletion-multiplicity exactly 2.

| conj | residual size | multiplicity of the worst vertex |
|---|---|---|
| 194 | 74 | `{1: 1, 2: 73}` |
| 198a | 30 | `{2: 30}` |
| 200 | 2 | `{1: 1, 2: 1}` |
| 217 | 56 | `{1: 5, 2: 51}` |

So each hypothesis, whatever else it does, silently forbids the one structure
that would refute it immediately. Proving that implication directly — for even
one of the four — would be a real step, since it removes the only obstruction
that is visible without going into the structure of the pieces.

Confirmed over all **12,111 connected graphs on 3–8 vertices**: 627 of them have
such a vertex, and **none satisfies any of the four hypotheses** (`_cutcheck.py`).

### Where a counterexample would have to be, and how close it gets

A counterexample to any of the four is a graph with **no** Hamiltonian path that
satisfies the hypothesis. Of the 994 connected graphs on 3–7 vertices, 144 have
no Hamiltonian path. Measuring how close those come:

| conj | best slack among non-Hamiltonian graphs | needs |
|---|---|---|
| 198a | `−0.1429` | `≥ 0` |
| 194 | `−0.2857` | `≥ 0` |
| 217 | `−1.0` | `≥ 0` |
| 200 | closest `\|tree − ⌈1+l_avg⌉\| = 1` | `= 0` |

`#198a` misses by one seventh, and these slacks move in steps of `1/n`, so the
granularity gets finer as `n` grows. That makes a targeted search worthwhile —
over graphs constructed to have no Hamiltonian path, maximising the hypothesis
slack — rather than over graphs in general.

### Applying the 322 lesson to the other conditionals

Measuring what each hypothesis admits among the 995 connected graphs on ≤ 7
vertices: `#322` admits 6 (and they are `K₂`–`K₇`, which is the whole content of
that conjecture), `#316` admits 27, `#314` admits 40. Neither of the latter two
collapses — they are narrow but genuinely varied, mostly bipartite for `#314`.

`#314` almost lands inside a known theorem. Bahadır, Ekim and Gözüpek
(*Well-Totally-Dominated Graphs*, [arXiv:2010.02341](https://arxiv.org/abs/2010.02341))
characterise **triangle-free WTD graphs with `γ_t = 2`**. Of the 40 graphs
satisfying `#314`'s hypothesis, **35 have `γ_t = 2`** and so fall under that
characterisation; the remaining 5 have `γ_t = 3` and do not. Not a resolution,
but it locates the gap precisely.

## What went wrong (continued)

**13. The machine was saturated, and two orphaned processes were burning cores
unnoticed.** The parallel product scan was launched with 14 workers on a
16-core machine, which is correct if the machine is a batch node and wrong if
somebody is sitting at it. They were, and everything they did became slow.

Two further processes were found running that should not have been. One was a
3-cell scan for bound 46 started nearly four hours earlier: it had been launched
under a `timeout`, was manually backgrounded instead, produced **no output at
all** in that time, and had accumulated 13,951 seconds of CPU — on a question
that `family46.py` had since answered analytically. The other was a
`ham_search.py` started with a shell `&`, which put it outside the tracking that
would have let it be stopped; a later, properly launched copy ran *alongside* it
rather than replacing it.

Three separate mistakes with one root: **treating compute as free**. It is not
free when it is someone's own machine.

The rules that follow, and are now followed:

* **cap workers well below core count** — 4 of 16 here. The gain from 4 to 14 is
  a few hours on a job that is not urgent; the cost is an unusable machine.
* **never launch with a bare `&`** — anything not launched through the tracked
  mechanism cannot be stopped, and a "replacement" run silently doubles the load.
* **make long jobs print progress.** The four-hour silent one would have been
  caught within minutes of its first missing report.

## Where the results actually came from

Worth tabulating, because it is not where the effort went.

| result | how it was found |
|---|---|
| `γ_t(G □ K₂) = γ(G × K₂)` | reading the **equality cases** of a scan |
| conjecture 322 resolved | measuring what the **hypothesis admits** |
| conjecture 100's prose/statement mismatch | **reading the definitions** |
| bounds 44 and 46 asymptotic families | reading **what the search converged to** |
| conjecture 160 false as formalised | exhaustive search — and it found a *transcription bug*, not mathematics |

Counterexample search produced one finding all day, and that one was a bug in a
Lean file. Everything else came from reading the structure of computations that
found nothing.

This is not surprising in hindsight. Searching for counterexamples puts you in
direct competition with groups running deep reinforcement learning and Monte
Carlo tree search over the same conjectures. Reading equality cases, hypothesis
strength and definitions competes with nobody, because it is not what those
methods produce.

## What went wrong (continued)

**11. Formalisation removes the ambiguity, not the misreading.** The whole point
of using the Lean list was that Graffiti's definitions are unreliable. The first
run of conjecture 144 then reported **25 counterexamples**, the first being the
triangle `K₃`.

The bug was mine: `ecc(S)` in the formalisation is the *set* eccentricity — the
maximum over vertices **outside** `S` of the distance **to** `S` — and it had
been implemented as "the eccentricity of the vertices in the centre", which is
just the radius. On `K₃` every vertex is a centre, so the true value is 0 and the
conjecture holds with equality; the wrong value made it fail.

It was caught only because the counterexample was a triangle. A three-vertex
counterexample to a conjecture on a published list is not a discovery, it is a
bug, and the size of the witness is the fastest thing to check. Formalisation
guarantees the statement is unambiguous. It guarantees nothing about whether the
reader implemented that statement.

## Every conjecture in range but one

The free-bound pass finished. Of the thirteen conjectures the toolkit tracks,
**twelve are now settled over all 995 connected graphs on at most 7 vertices by
bounds alone** — no conjecture-specific computation is needed to see that they
hold there.

```
 conj   covered  residual
    2       995         0
   19       995         0
   40       995         0
   59       995         0
   61       995         0
  133       995         0
  141       995         0
  142       995         0
  144       995         0
  145       995         0
  146       995         0
  291       994         0   (of 994 where the hypothesis applies)
  100       639       100   (of 739 where the hypothesis applies)
```

### The bounds added in this pass

All verified with zero violations over the 995 graphs.

| bound | reason |
|---|---|
| `tree(G) ≥ n − 1` when deleting one vertex leaves a tree | the deletion *is* the witness; the test is one connectivity-and-acyclicity check per vertex |
| `tree(G) ≥ 2·rad − 1` | a shortest path realising the radius, extended both ways |
| `tree(G) ≥ 3 + max_v l(v)` when `girth ≥ 6` | the star at `v` extends another step without closing a cycle |
| `f(G) ≥ n − (m − n + 1)` | delete one vertex per independent cycle |
| `Ls(G) ≥ Δ` | root a spanning tree at a maximum-degree vertex |
| `α ≤ n − ν` | Gallai: each edge of a maximum matching keeps a vertex out |
| `γ_t ≤ n − Δ + 1` | classical |
| `γ_t ≤ 2γ` | each dominator drags in one neighbour |

The cyclomatic bound alone took **conjecture 59 from 22 residuals to 2**.

Three of these run the *other* way — they bound the left-hand side from above —
so the toolkit now carries a mirrored table: prove `lhs ≤ rhs` by squeezing an
upper bound on `lhs` through. Conjectures 100 and 291 live there.

### Two conjectures that fell to a reading rather than a bound

**145** (`2·eccSet(B) ≤ tree(G)·lMin(Ḡ)`) splits cleanly on the one factor
nobody looks at. If `lMin(Ḡ) ≥ 2` the right-hand side is at least `2·tree`,
which is at least `2(diam+1)`, while `eccSet(B) ≤ diam` always — so it holds
**outright, for every graph**, with no computation. If `lMin(Ḡ) = 1` the
statement collapses to `2·eccSet(B) ≤ tree(G)`, which the existing tree
toolkit answers. Residual 0.

**291** (`γ_t ≤ HH-zero-step + freqMinTriangles`) closes on an exact witness
rather than a bound: `γ_t = 2` **exactly when some edge's open neighbourhoods
cover `V`**, which is an `O(mn)` check instead of a set-cover solve. Both
residual graphs had `γ_t = 2` against a right-hand side of 3; the general
bounds `n − Δ + 1` and `2γ` both overshot to 4.

### Witness certificates, and why they are a different kind of object

The last three residuals — conjecture 2's nine, 59's two, 133's one — fell to
objects that must be labelled honestly:

* `Ls ≥ n − |S|` for any connected dominating set `S` **greedy can exhibit**;
* `path(G) ≥ |P|` for any induced path **greedy can exhibit**;
* `f(G) ≥ |I ∪ A|` for `I` a maximum independent set and `A` independent with
  the bipartite graph between them **acyclic** — a genuine generalisation of the
  star-forest bound, whose condition (one neighbour in `I` each) is sufficient
  for acyclicity but not necessary.

A star, a ball, Gallai's identity — these are ingredients of a proof for *all*
graphs. A witness exhibits one object on the graph in hand and implies nothing
about any other graph. **They finish the in-range verification; they are not
steps toward a theorem**, and the file says so at the point of definition so the
distinction cannot quietly erode.

### The one that does not move, and why

**Conjecture 100 is stuck for a structural reason, not for want of trying.** It
bounds `α` from *above*. Every technique in this project works by exhibiting a
witness — an induced star, a ball, a greedy path, a connected dominating set —
and a witness proves a maximum-type invariant is **large**, never that it is
small. Bounding `α` from above means ruling out *all* independent sets, which is
computing `α`.

Four separate valid upper bounds were tried (`n − ν`, `n − ⌈(n−1)/Δ⌉`, `n/2`
under a perfect matching, `n − 1`) and the residual stayed at exactly 100 — the
extra three add nothing beyond Gallai. The residuals hold with **equality**,
which is the signature of a bound that is not merely weak but aimed wrong.

This is the sharpest thing the pass produced about its own method: **the
free-bound method is one-directional.** It settles lower bounds on hard
invariants cheaply and has nothing at all to say about upper bounds on them.
Nine conjectures fell to it; the one that points the other way did not budge.

## Result: conjecture 314 is true

> **WOWII 314.** *A connected triangle-free graph whose largest induced path has
> at most 4 vertices is well totally dominated — all its minimal total
> dominating sets have the same size.*

Previously untouched by this project. **Resolved.**

Throughout, `S` is a minimal total dominating set and `w(x)` is a *private
total-neighbour* of `x ∈ S`: minimality means `S − {x}` fails, so some `w(x)`
has `N(w(x)) ∩ S = {x}`. The `w(x)` are distinct, since that equation recovers
`x` from `w(x)`. **Every step is an induced `P₅`, produced from the structure
and killed by the hypothesis.**

### Part 1: `|S| ≤ 3`

Suppose `|S| ≥ 4`.

**(a) `G[S]` has no isolated vertex** — total domination asks every vertex of
`V` to have a neighbour in `S`, including the vertices of `S`.

**(b) `G[S]` is connected.** Let `C₁, C₂` be components, `k = dist(C₁, C₂)`.
There is no edge between them (`G[S]` is *induced*), so `k ≥ 2`; and `k ≥ 4`
would make a shortest path between them an induced `P₅` outright. Pick `x₁, y₁`
adjacent in `C₁` and `x₂, y₂` adjacent in `C₂`, which (a) permits.

* `k = 2`, with `z` a common neighbour of `x₁, x₂`: `z ∉ S`, since an `S`-vertex
  adjacent to both would join the components. Then `y₁ — x₁ — z — x₂ — y₂` is
  induced: `y₁z` and `y₂z` would close triangles, and `y₁x₂`, `y₁y₂`, `x₁y₂`,
  `x₁x₂` all lie across two components.
* `k = 3`, along a shortest `x₁ — u₁ — u₂ — x₂`: then `y₁ — x₁ — u₁ — u₂ — x₂`
  is induced — `y₁u₁` closes a triangle, `y₁u₂` would put `y₁` within distance 2
  of `x₂`, and `x₁u₂`, `u₁x₂` are excluded because the path is shortest.

**(c) `G[S]` has no induced `P₄`.** If `p₁—p₂—p₃—p₄` is induced in `G[S]`, then
`w(p₁) — p₁ — p₂ — p₃ — p₄` is induced in `G`: `p₂, p₃, p₄ ∈ S` and `p₁` is the
only `S`-neighbour of `w(p₁)`, which also forces `w(p₁)` to differ from each.

**(d) So `G[S] = K_{s,t}`.** A `P₄`-free graph is a cograph; a connected cograph
is a join; and triangle-freeness forces both sides of the join to be edgeless.

**(e) Neither side has three vertices.** With `a₁, a₂, a₃` on one side and `b`
on the other, the three privates cannot be pairwise adjacent — that is a
triangle — so some pair `w(a₁), w(a₂)` is not, and
`w(a₁) — a₁ — b — a₂ — w(a₂)` is induced.

**(f) So `s = t = 2` and `G[S]` is a 4-cycle `a₁—b₁—a₂—b₂—a₁`.** No `w` lies in
`S`, since every vertex of a 4-cycle has two `S`-neighbours while a private has
one — so the eight vertices are distinct. And now **every one of the 28 pairs is
forced**:

* `w(a₁)w(a₂)` is an edge, else `w(a₁) — a₁ — b₁ — a₂ — w(a₂)` is induced;
  likewise `w(b₁)w(b₂)`.
* `w(a₁)` is adjacent to `w(b₁)` or `w(b₂)`, else
  `w(a₁) — a₁ — b₁ — w(b₁) — w(b₂)` is induced. Say `w(a₁)w(b₁)`.
* Then `w(a₂)w(b₁)` is a non-edge (triangle with `w(a₁)`), so the previous point
  forces `w(a₂)w(b₂)`; and `w(a₁)w(b₂)` is a non-edge (triangle with `w(b₁)`).
* Everything else is a non-edge: an `S`-vertex meets no private but its own, and
  `a₁a₂`, `b₁b₂` are the diagonals.

That names **one graph**: 3-regular on 8 vertices with girth 4 — the **Wagner
graph `V₈`**, the circulant `C₈(1,4)`. And `V₈` contains an induced `P₅`, for
instance `b₁ — a₁ — b₂ — w(b₂) — w(a₂)`. Contradiction. So `|S| ≤ 3`. ∎

### Part 2: all minimal sets have equal size

A total dominating set has at least 2 vertices, and one of size 2 is exactly a
**dominating edge** `uv` with `N(u) ∪ N(v) = V`.

**If `G` has a dominating edge**, every minimal `S` has `|S| = 2`.
Triangle-freeness makes `N(u) ∩ N(v)` empty, so `A := N(u)` and `B := N(v)`
**partition** `V` and each is independent: `G` is bipartite, `u` adjacent to all
of `A`, `v` to all of `B`. Any total dominating set meets both parts. If
`a₁, a₂ ∈ S ∩ A` were distinct, their privates `w(a₁), w(a₂)` lie in `B` and

    w(a₁) — a₁ — u — a₂ — w(a₂)

is induced — `w(a₁)w(a₂)`, `w(a₁)u`, `uw(a₂)` are within `B` and `a₁a₂` within
`A`; `w(a₁)a₂` and `w(a₂)a₁` are excluded by privateness; and `u ≠ w(a₁), w(a₂)`
because `N(u)` is all of `A`. So `|S ∩ A| ≤ 1`, and running it with `v` gives
`|S ∩ B| ≤ 1`.

**If `G` has no dominating edge**, no minimal `S` has size 2, and Part 1 caps it
at 3. So `|S| = 3`. ∎

### Notes

**Bacsó–Tuza is not enough.** Their theorem — every connected `P₅`-free graph
has a dominating clique or a dominating `P₃` — gives `γ_t ≤ 3` across the class
immediately, since triangle-freeness shrinks the clique to a vertex or an edge
and each is already a total dominating set. That explains the *value* 3. It says
nothing about minimal sets overshooting, which is the entire content.

**Triangle-freeness is load-bearing.** Of the 637 connected `P₅`-free graphs on
at most 7 vertices, **209 fail** to be well totally dominated, the smallest on 5
vertices — edges `04, 12, 13, 23, 34`, with minimal sets of sizes 2 and 4. Every
use of the hypothesis above is a triangle closed by one edge.

**How it was found, and the one honest use of a computer.** By writing down what
a minimal set's private neighbours force and reading off the induced `P₅` each
time. One branch — `s = t = 2` — resisted; but it forced *every* adjacency among
eight vertices, so it named a single graph, and asking whether that graph
satisfies the hypothesis took a second. **The search space had already collapsed
to a point; the computation only confirmed it.** That is the opposite of the
counterexample searches that consumed most of this project's compute and
produced one bug.

## Result: conjecture 100 holds whenever `α ≥ 16`

> **WOWII 100.** *`α(G) ≤ ⌈(max_v l(v) + 0.5·‖deg(Ḡ)‖₂)/2⌉`.*

The one conjecture the free-bound method could not touch, because it bounds a
hard invariant from *above*. It yields to attacking the **other side** of the
inequality instead.

> **Lemma 1.** `‖deg(Ḡ)‖₂ ≥ (α−1)√α`.
>
> *Proof.* A maximum independent set `I` of `G` is a **clique of size `α` in
> `Ḡ`**, so each of its `α` vertices has `Ḡ`-degree at least `α−1`. Hence
> `Σ_v deg_Ḡ(v)² ≥ α(α−1)²`. ∎

> **Lemma 2.** `max_v l(v) ≥ 2` for every connected non-complete graph.
>
> *Proof.* Some pair is non-adjacent, so a shortest path between such a pair has
> at least three vertices; its first three give `x, z, y` with `x, y ∈ N(z)` and
> `xy` a non-edge, so `l(z) ≥ 2`. ∎

Since `⌈t⌉ ≥ t`, the right-hand side is at least `1 + (α−1)√α/4`, and

    α ≤ 1 + (α−1)√α/4   ⟺   (α−1) ≤ (α−1)√α/4   ⟺   4 ≤ √α   ⟺   α ≥ 16.

**Theorem.** Conjecture 100 holds for every graph with `α(G) ≥ 16`; it can only
fail in the range `α ≤ 15`. ∎ (The threshold is exact for this argument —
equality at `α = 16`.)

**Why the free bounds failed and this did not.** Every other technique here
exhibits a **witness**, and a witness proves a maximum-type invariant is
*large*, never that it is small. Four valid free upper bounds on `α` were tried
— Gallai's `n − ν`, `n − δ`, the clique bound `(1+√(1+8m̄))/2`, and `n/2` under a
perfect matching — and together they left **98** of the 739 applicable graphs
untouched, against **100** for Gallai alone. The other three add almost nothing.

What works is not bounding `α` at all but showing the **right-hand side grows
with `α`** — which it must, because a large independent set is a large clique in
the complement and therefore drives complement degrees up. Same move as
everywhere else in this project, aimed at the other side.

**What is left.** For `α ≤ 15` the proof throws away a real trade-off. Let `t`
be the largest number of `I`-neighbours of a vertex outside `I`. Then
`max_v l(v) ≥ t`, since `N(u) ∩ I` is independent; and every vertex outside `I`
has `Ḡ`-degree at least `α − t`, so

    ‖deg(Ḡ)‖₂² ≥ α(α−1)² + (n−α)(α−t)².

Large `t` helps the first term of the right-hand side, small `t` helps the
second. Working that seesaw is what would push the threshold below 16.

## The four Hamiltonian-path conjectures — and 200 is false

> **194.** `α ≤ 1 + l_avg` ⟹ Hamiltonian path
> **198.** `b(G) ≤ 2 + avg ecc` ⟹ Hamiltonian path
> **200.** `tree(G) = ⌈1 + l_avg⌉` ⟹ Hamiltonian path
> **217.** `Ls(G) ≤ 4·[residue = 2] + 2` ⟹ Hamiltonian path

All four have the same shape: an invariant is small, therefore a Hamiltonian
path exists. So the question is not *can we find a Hamiltonian path* but **which
classical sufficient condition does the hypothesis already imply**.

Five certificates, each a theorem and each checked against `has_hamiltonian_path`
over all 995 graphs before being used — a subtly misimplemented certificate
looks exactly like a strong one:

| certificate | note |
|---|---|
| Chvátal–Erdős `κ ≥ α − 1` | |
| Ore `deg(u)+deg(v) ≥ n−1` | |
| `Ls(G) ≤ 2` | every spanning tree has ≥2 leaves, so the maximum being 2 makes **every** spanning tree a path |
| **Bondy–Chvátal closure of `G ∨ K₁`** | a Hamiltonian path in `G` is a Hamiltonian cycle in `G ∨ K₁` |
| Chvátal's degree condition on `G ∨ K₁` | |

```
 conj   applicable   residual
  194          501          1
  198          148          3
  200          136          0   <-- settled in range
  217          453          0   <-- settled in range
```

The same run over the 12,112 connected graphs on at most **8** vertices:

```
 conj   applicable   residual
  194         5621         19
  198          516         16
  200          671          0   <-- still settled
  217         4013          0   <-- still settled
```

`200` and `217` hold up: 671 and 4,013 applicable graphs, and classical theory
covers every one. That is what makes the proof target credible rather than a
small-case coincidence.

**The closure is what does the work.** Ore, Dirac and Chvátal–Erdős are one-shot
degree or connectivity tests. The Bondy–Chvátal closure *iterates* — each edge it
adds raises degrees and lets it add more — and it alone covers 497 of 194's 501,
144 of 198's 148, and all of 200's and 217's.

### What that says about the conjectures

**These hypotheses are not saying anything new about Hamiltonicity.** They are
roundabout ways of asserting a condition Bondy and Chvátal covered in 1976.

> ### ⚠ The proof target proposed here for 200 cannot exist
>
> This section originally concluded that the thing to prove was
> `hypothesis ⟹ the closure of G ∨ K₁ is complete`, and that for **200** the
> sharper target `tree(G) = ⌈1 + l_avg⌉ ⟹ Chvátal's condition on G ∨ K₁` would
> do it. **Conjecture 200 is false.** Jitendra Prajapati found an 11-vertex
> counterexample on 21 July 2026, `graph6 J??FFBRq}N_`, verified in
> `wowii_status.py`: `tree(G) = ⌈1 + l_avg⌉ = 4` and there is no Hamiltonian
> path. The exhaustive range here was 8 vertices, so nothing in the computation
> could have seen it.
>
> The certificates themselves were **not** fooled — on that graph all five
> correctly decline. The tooling was sound; the literature check was not.

For **217**, still open, the closure target stands.

**The easy obstruction never appears.** No graph satisfying any of the four
hypotheses has a vertex cut of size `k` leaving more than `k+1` components — the
standard reason a graph has no Hamiltonian path — checked for `k ≤ 3` over every
graph in range. The hypotheses rule it out on their own, which is consistent
with all four being true.

The 4 residual graphs are all `n = 7`, `m ∈ {7,8}`, `κ = 1`: sparse, nearly
unicyclic, one cut vertex. Degree-based certificates are blind to exactly this
shape, and a block-decomposition certificate is the natural next tool.

## Verified over all graphs on 8 vertices

Both conjecture batches were rerun over the **12,112 connected graphs on 2..8
vertices**, generated by `allgraphs.py`, which count-checks itself against OEIS
A000088 and A001349 before returning anything.

```
#2    Ls(G) >= 2*(l_avg - 1)                        applicable 12112  holds 12112
#19   floor(avg ecc + max l(v)) <= b(G)             applicable 12112  holds 12112
#40   f(G) >= ceil((p(G)+b(G)+1)/2)                 applicable 12112  holds 12112
#59   f(G) >= ceil(sqrt(residue*b))                 applicable 12112  holds 12112
#61   f(G) >= residue(G) + ceil(diam/3)             applicable 12112  holds 12112
#65   distMin(A) + ceil(distMin(M)/3) <= f(G)       applicable 12112  holds 12112
#100  alpha <= ceil((max l(v) + 0.5*||deg Gc||)/2)  applicable 10627  holds 10627
#141  tree(G) >= floor(girth/2) - 1 + max l(v)      applicable 12112  holds 12112
#144  tree(G) >= girth - 1 + ecc(center)            applicable 12112  holds 12112
#145  2*eccSet(B) <= tree(G) * lMin(complement)     applicable 10860  holds 10860
#194  alpha <= 1 + l_avg  =>  Hamiltonian path      applicable  5621  holds  5621
#198  b(G) <= 2 + avg ecc  =>  Hamiltonian path     applicable   516  holds   516
#200  tree(G) == ceil(1+l_avg)  =>  Hamiltonian path applicable  671  holds   671
#291  gamma_t <= HH-zero-step + freq(min triangles) applicable 12111  holds 12111
#314  triangle-free, induced path <= 4  =>  WTD     applicable    81  holds    81
#316  avgdeg(complement) <= #pendants  =>  WTD      applicable    37  holds    37
#322  n>=5, every l(v)<=1  =>  WTD                  applicable     4  holds     4
```

No counterexample at `n = 8`. This is verification, not evidence of much — the
value is that it would have caught a transcription error, which is exactly how
conjecture 160's Lean bug surfaced.

## What went wrong (continued)

**14. A tidy-up dropped a bound because one looked more general.** Consolidating
the scattered bounds into `wowii_toolkit.py` replaced the single-path cycle bound
with the caterpillar version, on the assumption that the caterpillar generalises
it. **It does not, in either direction.** The caterpillar attaches single
vertices at *several* cycle vertices; the path version attaches a *whole path* at
*one*. Neither contains the other.

Conjecture 144 regained a residual and stayed there until both were present, and
its one hard case — a `C₅` with a length-2 path hanging off a single vertex — is
visible only to the path version. "The more general-looking one subsumes the
other" is a guess, and it was made silently while cleaning rather than while
proving.

**15. A batch test was contaminated by unvalidated bounds.** A candidate sweep
reported conjectures 133, 142 and 146 all dropping to residual 0. Two of those
were real. The third was not: the sweep applied *every* candidate to the
residuals, including two that had **just been printed as invalid** in the same
output — one with 567 violations, one with 182.

The validity column and the coverage column were computed in the same script and
one did not gate the other. The fix is not to be more careful reading output; it
is that a bound must be *rejected* before the coverage loop can see it. Caught
within one step by rerunning the real toolkit, whose bound list only ever
contains validated entries — which is an argument for the toolkit being the only
path to a coverage number.

## The literature check, done properly — and what it found

The Lean formalisations carry `category research open`. **That tag is not a
status report.** It records what was believed when the file was written. The
authoritative source is DeLaViña's own page, split into `resolved.htm` and
`open.html`, and **still being updated** — two of these conjectures were refuted
in the week of 21–23 July 2026, days before this check.

Checked 27 July 2026. The page's numbering was matched to the Lean numbering by
comparing the **statement text of all 22**, not the numbers.

| conj | status per the primary source |
|---|---|
| **146** | **PROVED**, 21 Jul 2026, Brain Akaka — a human-readable proof, developed and checked with AI systems, upstreamed as PR 4505 to `formal-conjectures` |
| **198** | **RESOLVED**, June 2010, Richard Stong — but see the note below: the statement implemented here is **198a**, a *different* conjecture, still open |
| **200** | **FALSE**, 21 Jul 2026, Jitendra Prajapati — 11-vertex counterexample, `J??FFBRq}N_` |
| **291** | **FALSE**, 23 Jul 2026, Zyad Tamimi — 12-vertex counterexample |
| the other 18 | still listed as open |

Both refutations are verified in `wowii_status.py`. **291 was reproduced
independently**: hill-climbing on `γ_t − (k + freq)` found a *different*
12-vertex graph with the same profile (`γ_t = 4, k = 2, freq = 1`) and none on
11 vertices, consistent with 12 being the minimum order.

### What it cost

**198 versus 198a — a reprieve, and a warning.** The page carries both.
Conjecture **198** asks `b(G) ≤ 2 + ecc_avg(M)` with `M` the *maximum-degree*
vertices; Stong resolved it in 2010. Conjecture **198a** asks
`b(G) ≤ 2 + ecc_avg(G)`, averaged over *all* vertices, and is **still open**.
`conj_198a` implements the second, so that effort was not spent on a settled
problem. But it was recorded under the label "198" everywhere — in the code, the
tables, the commit messages — which is exactly how such a confusion propagates,
and it took reading the source page to notice that the label and the code
disagreed.

**DeLaViña's own notes carry partial results, and one of them is a bound this
project derived independently.** Under conjecture 40, dated 6 March 2004:
`f(G) ≥ b(G)/2 + 1`, with the corollary that 40 follows whenever the path
covering number is 1 — which is 851 of the 995 graphs in range. It is valid
(zero violations) and is now in the toolkit, attributed. Under conjecture 19,
dated 22 June 2005: if `ecc_avg ≤ diam − 1` then 19 follows from WOWII 13 —
though that covers only 3% of the graphs in range, well short of what the free
bounds do.

The lesson compounds: **the status page was not the only thing left unread. The
conjecture entries themselves carry twenty years of accumulated notes**, and
they are the first place a partial result would already be recorded.

Effort went into **146**, proved six
days earlier — by someone using the same combination this project uses, a
human-readable argument developed with AI assistance. And into pushing **200**
and **217** as a pair, concluding that the proof to attempt was
`hypothesis ⟹ the Bondy–Chvátal closure is complete`. For 200 that proof cannot
exist.

The consolation is narrow but real: **the tooling never lied.** On Prajapati's
graph all five Hamiltonicity certificates correctly decline, and the in-range
claims (`n ≤ 8`) were all literally true. The failure was entirely in not
reading the source list.

## Result: conjecture 142 is true for girth at most 3, by the same move

> **WOWII 142.** *`tree(G) ≥ (2/3)·girth(G) + ecc(B)`, with `B` the periphery —
> the vertices of **maximum** eccentricity.*

> **Lemma.** `diam(G) ≥ 1 + ecc(B)`.
>
> *Proof.* Let `v` attain `k = ecc(B) = dist(v, B)`. If `k = 0` the claim is
> `diam ≥ 1`. If `k ≥ 1` then `v ∉ B`, and `B` is exactly the set of vertices
> whose eccentricity **is** the diameter, so `ecc(v) ≤ diam − 1`; and some
> vertex of `B` sits at distance `k`, so `k ≤ ecc(v) ≤ diam − 1`. ∎

One line, where the centre version needed two — the periphery is defined by the
diameter, so both halves of the argument collapse together. Then for a graph
containing a triangle,

    tree(G) ≥ diam + 1 ≥ 2 + ecc(B) = (2/3)·girth + ecc(B),

and forests are trivial. **142 holds for every connected graph of girth ≤ 3** —
11,803 of the 12,112 graphs on at most 8 vertices.

**The two lemmas are worth stating on their own.** For any connected graph on
`≥ 2` vertices, `diam ≥ 1 + ecc(Centers)` and `diam ≥ 1 + ecc(Periphery)`. Both
have zero violations over the 12,112-graph corpus, neither appears in the
literature reached here, and each buys a whole girth regime of a conjecture that
had no proof at all.

## Result: conjecture 144 is true for girth at most 3

> **WOWII 144.** *`tree(G) ≥ girth(G) − 1 + ecc(Centers(G))`.*

The whole regime turns on one elementary fact about centres that appears not to
be written down anywhere.

> **Lemma.** `diam(G) ≥ 1 + ecc(Centers(G))` for every connected graph on at
> least two vertices.
>
> *Proof.* Let `C` be the centre, `r = rad(G)`, and let `v` attain
> `k = ecc(C) = dist(v, C)`.
>
> If `k = 0` the claim is `diam ≥ 1`, true since `G` is connected on `≥ 2`
> vertices.
>
> If `k ≥ 1` then `v` is **not** a centre, so `ecc(v) > r` and hence
> `diam ≥ ecc(v) ≥ r + 1`. And for any `c ∈ C`, `dist(v,c) ≤ ecc(c) = r`, so
> `k ≤ r`. Therefore `diam ≥ r + 1 ≥ k + 1`. ∎

**Corollary.** If `G` contains a triangle then

    tree(G) ≥ diam + 1 ≥ 2 + ecc(Centers) = girth − 1 + ecc(Centers),

the first step because a shortest path is an induced tree. Forests are trivial
(`tree = n`, right-hand side `ecc(C) − 1`). So **144 holds for every connected
graph of girth ≤ 3** — 930 of the 995 graphs in range, and a complete proof for
the regime rather than a verification.

Zero violations of the lemma over the 12,112 connected graphs on ≤ 8 vertices.

**What is left, precisely.** Girth `≥ 4`, where the target `(g−1) + ecc(C)`
outruns what a shortest path supplies. The obstruction is concrete: **a shortest
cycle need not be eccentric enough.** On 24 of the 309 graphs of girth ≥ 4 on at
most 8 vertices — the smallest on 6 — *no* shortest cycle `C` has
`ecc(C) ≥ ecc(Centers)`, so the cycle-path construction cannot reach the target
on its own. The conjecture measures eccentricity from the **centre**; every
construction here starts at the **cycle**. Closing it means relating those two,
or building from the centre outward.

## How likely is any of the rest to be proved

Asked directly, and worth answering with numbers rather than encouragement.

### The base rate

DeLaViña's page splits into `resolved.htm` (**250** conjectures) and
`open.html` (**160**). So roughly **61% have fallen in twenty-odd years**, and
what is left is the residue that survived that filtering. A conjecture still
open in 2026 is open because the easy things were tried.

Against that, the list is **not dormant**. Two conjectures were refuted in the
week of 21–23 July 2026 and one was proved on 21 July — and that proof is
described on the page as "developed and checked collaboratively with AI
systems". Whatever this project is doing, other people are doing it too, right
now, and reaching the same conjectures within days.

### What the coverage data actually says

The honest split is **not** "settled in range" versus "open". It is whether the
*structural* bounds — the ones that are proof material for all graphs — keep up
as `n` grows. Rerunning the toolkit with the greedy witnesses removed:

| conj | structural residual, `n ≤ 7` | structural residual, `n ≤ 8` |
|---|---|---|
| 19, 40, 61, 141, 142, 144 | 0 | **0** |
| 133 | 1 | 1 |
| 59 | 2 | 14 |
| 2 | 9 | 32 |
| 145 | 0 | **75** |
| 146 | 0 | **96** |

That table is the whole answer.

**Where the residual stays 0, a proof is plausible.** The bounds already cover
every graph in reach; what is missing is an argument that the case analysis
exhausts *all* graphs, which is exactly the shape of the proof that worked for
141 — two regimes, split at girth 6, and the split point forced. Six
conjectures sit here. This is the range where I would expect real progress, and
where I would put maybe **one in three** on a complete proof for any given one
with sustained effort.

**Where the residual grows, the bounds are not the proof.** For 145 and 146 the
coverage was 0 at `n ≤ 7` and 75 and 96 at `n ≤ 8` — the free bounds simply
stop tracking the conjecture. Reading that as "nearly proved" would have been a
mistake, and it was the one thing the earlier optimism here got wrong. These
need a genuinely different idea, and **146 has now been proved by someone else
anyway**, which is the most direct evidence available that the free-bound route
was not the one.

**The Hamiltonian family should be attacked as false, not proved.** 200 held for
every one of the 671 applicable graphs on at most 8 vertices, was covered
completely by classical Hamiltonicity theory in that range, and is **false at
11**. That is a specific warning about 194, 198a and 217: the same hypotheses,
the same family, the same behaviour in the same range. `wowii_hunt.py` searches
9–14 vertices for exactly this.

### The uncomfortable part

**"Settled over all graphs on at most 8 vertices" carries very little
information about truth.** 200 is the proof: 12,112 graphs, every applicable
one holding, five independent classical certificates covering all of them — and
the statement is false. 291 likewise, at 12 vertices. The exhaustive range here
is not merely small, it is *below the scale at which these conjectures decide
themselves*, and the counterexamples that exist were found by targeted search,
not enumeration.

So the realistic breakdown of the 22:

* **6 moved here** — 65, 141, 314, 316, 322 proved; 160 refuted as formalised.
  Of these, 65 and 322 are trivial once read, 316 is a short counting argument,
  141 is real but its ingredients are textbook, and 314's `γ_t = 2` half turns
  out to be published. **The genuinely new mathematics is Part 1 of 314** — the
  cap `|S| ≤ 3` — and the `α ≥ 16` reduction of 100.
* **4 moved by others** — 146, 198, 200, 291.
* **12 left.** Of those, six look reachable by finishing a case analysis, three
  are Hamiltonian-family statements that may well be false, and three (2, 59,
  145) need something new.

### What I would do with the next hundred hours

Not more exhaustive computation — the range cannot be pushed far enough to
matter, and 200 shows what that range is worth. Two things instead:

1. **Finish one case analysis properly**, on 144 or 142, where the structural
   bounds already cover everything in reach. One complete proof is worth more
   than twelve partial ones.
2. **Hunt counterexamples above `n = 8` for the Hamiltonian family.** A
   counterexample is a definitive result, it is cheap when it exists, and the
   base rate in that family is now 1 for 1.

The thing not to do is what consumed most of this project's compute: exhaustive
scans over ranges too small to decide anything, and free bounds accumulated past
the point where they were still saying something.

## Attribution of every bound in the toolkit

Run at the same time as the status check, and it moved several entries out of
the "found here" column. **None of the bounds below were invented here except
where the last column says so**, and two that had been treated as fresh
observations are theorems from 1986.

| bound | due to |
|---|---|
| `path(G) ≥ 2·rad − 1` | **Erdős–Saks–Sós 1986, Thm 2.2** — proof by Fan Chung |
| `tree(G) ≥ 2·rad − 1` | same, since an induced path is an induced tree |
| `f(G) ≥ α + 1`, `tree ≥ 1 + max_v l(v)` | folklore — the neighbourhood star |
| `f, tree, path ≥ diam + 1` | folklore — a shortest path is induced |
| `tree(G) ≥ girth − 1` | folklore — a shortest cycle is chordless |
| `f(G) ≥ n − (m − n + 1)` | folklore; the systematic study of `t(G)` against the cyclomatic number is Erdős–Saks–Sós §3–7 |
| `b(G) ≥ f(G) ≥ tree(G) ≥ path(G)` | definitional |
| `Ls(G) ≥ Δ` | classical |
| `α ≤ n − ν` | **Gallai** |
| `γ_t ≤ n − Δ + 1`, `γ_t ≤ 2γ` | classical |
| Chvátal–Erdős, Ore, Bondy–Chvátal closure, Chvátal's degree condition | classical |
| Bacsó–Tuza (dominating clique or `P₃` in a `P₅`-free graph) | **Bacsó–Tuza 1990** |
| a dominating edge in a triangle-free graph splits `V` into two independent sets | **Bahadır–Ekim–Gözüpek 2020, Lemma 3.7** |
| `tree(G) ≥ Δ + ⌊girth/2⌋` for `girth ≥ 6` | assembled here from the standard locally-tree-like ball argument |
| `tree(G) ≥ 2 + max_v l(v)` for `girth ≥ 5` | assembled here |
| `f(G) ≥ |I ∪ A|`, `I` maximum independent, `A` independent, the bipartite graph between them acyclic | assembled here |

Erdős–Saks–Sós also prove (Thm 8.2) tight bounds for `t(G)` in terms of `α(G)`,
and their central result is that `t(G)` is **surprisingly small** for dense
graphs — `f(n, cn) = 2 log log n + O(log log log n)`. That is worth knowing
before spending effort on lower bounds for induced trees: there is very little
room in the dense range, and every bound here that works does so by living in
the sparse or high-girth range.

### What this means for conjecture 314

The dominating-edge half of the proof rests on a bipartition lemma that is
**Lemma 3.7 of Bahadır–Ekim–Gözüpek (2020)**, whose proof is the same three
lines. Worse for novelty, their **Theorem 3.8** characterises triangle-free
`WTD(2)` graphs completely, so that half of 314 is derivable from published
work — the route here is shorter for this particular class, but it is not new
mathematics.

**Part 1 — the cap `|S| ≤ 3` — does appear to be new.** Their paper handles
`γ_t = 2` in detail and poses "what are `WTD(k)` graphs" as an open problem; the
`γ_t = 3` case for triangle-free `P₅`-free graphs is not in it. That is the part
of 314 worth keeping: the private-neighbour argument that walks `G[S]` down to
complete bipartite and then names the Wagner graph.

## What went wrong (continued)

**16. Treating a formalisation's `open` tag as a status report.** The whole
reason for using the Lean list was that Graffiti's own statements are ambiguous,
and formalisation fixes that. It does — and it was allowed to answer a question
it was never claiming to answer. `category research open` says the statement was
open when someone wrote the file. It is not a live index.

The primary source was one HTTP request away and had a page literally named
`resolved.htm`. It says 198 was settled in **2010**. It carries dates from six
days before this check.

This is failure **5** again in a new costume. That one was "search for work
citing the conjecture list, not just the conjectures." The generalisation is
harsher: **for any curated list of open problems, find who maintains it and read
their status page first, before any mathematics.** Not at the end as a check —
first, as the cheapest possible filter. Twenty minutes of reading would have
redirected the whole Hamiltonian-path effort.

**17. A conjecture was proved in a weaker form than it was stated.** The
formalisation of 141 uses `⌊girth/2⌋` where DeLaViña's original has
`(1/2)·girth`. For odd girth the original is stronger, and the proof written
here did not reach it — on `C₇` it delivered 4 against a requirement of 5.

What makes this one instructive is that **no amount of computation could have
caught it.** Both readings were verified over all 12,112 connected graphs on at
most 8 vertices with zero violations, because there are exactly two graphs in
that range with odd girth ≥ 7, and the bound is slack on both. The two
statements are computationally indistinguishable in reach and mathematically
different. Only reading the original found it.

The gap turned out to be closable — see `wowii141.py` — so the result survives,
stronger than before. But it was proved against a transcription, not against the
conjecture.

## Files

```
search.py             cross-entropy search over edge indicators (kept, though it
                       lost to brute force on the validation target — should
                       matter at larger n)
exhaustive_trees.py   exhaustive enumeration over non-isomorphic trees, with a
                       linear-time tree matching routine
laplacian_bounds.py   all 68 conjectured bounds on μ, with the regular-graph
                       self-test that guards the transcription
known_bounds.py       classical proved bounds on μ, numerically re-validated
collatz_wielandt.py   the Lemma 3.6 machine behind most of the published proofs,
                       validated by reproducing them and by failing on every
                       refuted bound
dominance.py          searches for a classical theorem implying each open
                       conjecture pointwise, and the max-degree criterion
proofs.py             the twelve conjectures proved true and the inequalities
                       they rest on, re-verified as a regression test
family44.py           the near-extremal one-parameter family for bound 44 and
                       its 60-digit asymptotics
family46.py           the same for bound 46: semiregular bipartite, slack 1/(4k)
                       its 60-digit asymptotics
quotient_search.py    counterexample search over equitable partitions for
                       bounds 44 and 46, validated against the twelve published
                       counterexamples
sweep.py              parameter sweep driver for that search
zero_forcing.py       Z(G) and alpha(G), validated against ten published values
cubic_graphs.py       exhaustive generation of connected cubic graphs, with the
                       count check against OEIS A002851 that must pass first
scan_cubic.py         the conjecture over every connected cubic graph
allgraphs.py          exhaustive generation of all connected graphs, count-checked
                       against OEIS A000088 and A001349; lifts the 7-vertex cap
wowii.py              the untouched Graffiti.pc conjectures, from the Lean list
wowii322.py           conjecture 322 resolved, and why its hypothesis is vacuous
wowii_toolkit.py      every free bound, with its one-line reason, and the
                       residual each conjecture is left with
wowii141.py           conjecture 141 proved: neighbourhood star for girth <= 5,
                       locally-tree-like ball for girth >= 6
wowii61.py            conjecture 61 settled for diameter <= 3, via f >= alpha + 1
wowii316.py           conjecture 316 resolved by characterising its hypothesis
wowii160.py           the five-vertex counterexample to conjecture 160 as
                       formalised, and why no reading of c_C4 rescues it
domination.py         exact gamma and gamma_t by branch and bound, validated
scan_products.py      the product conjecture over all pairs of small graphs
search_products.py    local search over pairs, reaching asymmetric shapes the
                       exhaustive scan never sees
prism_identity.py     the H = K_2 theorem, with the bijection checked on subsets
cubic_search.py       local search over cubic graphs, for orders the exhaustive
                       scan cannot reach
triangle_replaced.py  a search-space reduction that turned out to be wrong;
                       kept as the record of failure 9
```

## Requirements

Python 3, `numpy`, `networkx`, `mpmath` (for the 60-digit asymptotics).
