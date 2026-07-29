# Counterexample search for graph-invariant inequalities

Status: amateur/independent investigation. Started 2026-07-26.

This repo attacks published open conjectures of the shape `f(G) ≥ g(G)` for
every (connected) graph `G` — a class where a counterexample is one finite
graph and verification is a calculation, following the approach of Wagner
(*Constructions in combinatorics via neural networks*, [arXiv:2104.14516](https://arxiv.org/abs/2104.14516)).

**The full research log — every attempt, every dead end, every retraction,
every "what went wrong" — is [RESEARCH_LOG.md](RESEARCH_LOG.md).** This file
is only the current bottom line: what's proved, what's refuted, what's still
open, and where to look in the code for each. Written in the order results
were found, the log contradicts itself in places as later sections correct
earlier ones — several headline claims below only became true after a
retraction. Treat this README, not the log, as the source of truth for
"where things stand now."

## 1. Validation target — resolved (reproduced, refined by one vertex)

**Conjecture** (Aouchiche et al. / Wagner's Conjecture 2.1): for connected `G`
on `n ≥ 3` vertices, `λ₁ + μ ≥ √(n−1) + 1`. Known false; Wagner found a
19-vertex counterexample.

**Result:** exhaustive enumeration of all non-isomorphic trees (`exhaustive_trees.py`)
finds the true minimum counterexample at **18 vertices**, one smaller than
reported, and proves no counterexample exists below it (adding an edge only
increases `λ₁ + μ`, so trees suffice). Used here mainly to shake out
tooling bugs before pointing anything at an open problem — see the log for
two false starts (a search that never found it, a floating-point ghost at
n = 14) that shaped every guard used afterward.

## 2. 68 Laplacian spectral-radius bounds — 12 resolved here, but scooped

**Target:** a family of 68 automatically-generated conjectured upper bounds
on the Laplacian spectral radius `μ(G)` (Brankov–Hansen–Stevanović). 36 were
still open when this project started.

**Result:** 12 of the 36 follow in one line each from classical theorems
(Anderson–Morley, Merris, Li–Zhang) via a max-degree argument — `dominance.py`,
`known_bounds.py`. **All twelve turned out to already be published**, six
weeks earlier, in Damnjanović–Ha–Stevanović ([arXiv:2606.14550](https://arxiv.org/abs/2606.14550)),
which leaves exactly **two conjectures open: #44 and #46**.

For those two: no counterexample found by exhaustive search (≤7 vertices,
≤18-vertex trees) or by a quotient-matrix search covering unbounded graph
size (`quotient_search.py`, validated against 12 published counterexamples).
Both bounds have an explicit family approaching equality and never crossing
it (`family44.py`, `family46.py`, to 60 digits) — evidence they are true, not
a proof. Three proof attempts, all recorded as failing, in the log.

**Status: #44 and #46 open. Everything else in this family is settled**
(22 confirmed true by this project or predecessors, 44 refuted by others).

## 3. Zero forcing vs. independence in cubic graphs — open, in progress

**Conjecture** (TxGraffiti / Davila et al.): connected cubic `G ≠ K₄` has
`Z(G) ≤ α(G) + 1`. No exhaustive check was previously reported.

**Result:** verified on every connected cubic graph up to 12 vertices (112
graphs); `K₄` is the only violation, as the conjecture says. `cubic_graphs.py`
generates the graphs (count-checked against OEIS A002851), `zero_forcing.py`
computes `Z`. Larger orders explored by local search (`cubic_search.py`), no
counterexample. **Open — search ongoing above n = 12.**

## 4. TxGraffiti product conjecture — proved for H = K₂

**Conjecture** (TxGraffiti): `γ_t(G □ K₂) ≥ γ(G × K₂)` for connected `G, H`.

**Theorem** (proved, `prism_identity.py`): equality holds for **every** graph
`G` when `H = K₂` — a size-preserving bijection between total dominating sets
of the prism and dominating sets of the bipartite double cover. This
generalises a 2017 result of Azarija–Henning–Klavžar (who needed `G`
bipartite) to every graph, with a shorter proof.

The general conjecture (arbitrary `H`) remains **open** — no counterexample
over 131,279+ pairs of small graphs, tight but never violated, and the `K₂`
proof technique is shown in the log not to extend past `H = K₂`.

## 5. 22 Graffiti.pc conjectures, from their Lean formalisation — the main result

**Source:** Google DeepMind's [`formal-conjectures`](https://github.com/google-deepmind/formal-conjectures),
`WrittenOnTheWallII` directory — 22 of DeLaViña's Graffiti.pc conjectures,
chosen because the Lean statement removes the definitional ambiguity that
sank earlier attempts on Graffiti-sourced conjectures.

Current status of all 22 (checked against DeLaViña's own `resolved.htm` /
`open.html`, not just the Lean `open` tag, which can be stale — two of these
were resolved by others in the week before this check):

| # | statement (short) | status |
|---|---|---|
| 65 | `distMin` bound | **open** — earlier "proved, trivial" used the wrong `distMin`; corrected 29 Jul 2026 after external review, no longer trivial, no counterexample found in range |
| 316 | avg-degree-of-complement bound ⟹ WTD | **proved** |
| 322 | `λmax(Ḡ) ≤ 1` ⟹ WTD | **proved** (complete multipartite) |
| 141 | girth/local-independence tree bound | **proved**, original (stronger) form |
| 314 | triangle-free + short induced path ⟹ WTD | **proved** |
| 100 | independence number upper bound | proved for `α ≥ 16`; open below |
| 19 | ecc/local-independence ⟹ bipartite number | proved in range (≤8 vertices); reduces to one open inequality |
| 40, 61 | forest-number bounds | proved in range; each reduces to 2 explicit bounds |
| 133 | radius/local-independence path bound | proved in range (995 graphs, n≤7) |
| 142, 144 | girth/eccentricity tree bounds | proved for girth ≤ 3 (144 also girth 4 mostly); small residual otherwise |
| 145 | eccentricity-set tree bound | proved except one explicit 6-vertex shape |
| 146 | eccentricity/radius-squared tree bound | **proved by someone else** (Brain Akaka, 21 Jul 2026) |
| 160 | 4-cycle/local-independence spanning-tree bound | Lean transcription is **false** (5-vertex counterexample); DeLaViña's actual definition holds except 11 residual graphs |
| 194, 198a, 217 | small-invariant ⟹ Hamiltonian path | mostly reduce to Chvátal–Erdős / Bondy–Chvátal closure; **217's open content is confined to n ≤ 10**; 194/198a small residuals |
| 198 | (same shape, different invariant than 198a) | **already resolved**, Stong 2010 — not attempted here |
| 200 | (same family) | **false** — 11-vertex counterexample (Jitendra Prajapati, 21 Jul 2026) |
| 291 | domination-number bound | **false** — 12-vertex counterexample (Zyad Tamimi, 23 Jul 2026); reproduced independently here |
| 2 | leaves/local-independence bound | proved for triangle-free graphs; reduces to one open inequality otherwise |
| 59 | forest/residue/bipartite bound | reduces to an AM–GM strengthening, proved in range |

**Net: 4 conjectures proved outright here (141, 314, 316, 322), several
more reduced to a single explicit inequality or a small finite residual, 2
resolved by others in the same weeks (146, 200-as-false), 1 already resolved
under a different label (198), and 291 independently confirmed false.** 65
was briefly counted as a fifth proof but that rested on a misreading of
`distMin`; see below. All the exhaustive claims (other than 65's, pending a
rerun) are verified over every connected graph up to 8 vertices (12,112 graphs,
count-checked against OEIS) via `wowii_toolkit.py` / `wowii_status.py`.

Two genuine Lean-transcription bugs were found and reported upstream
(conjectures 141, 322); a third (160) had already been fixed a day before
this check confirmed it independently. See the log for how each was found —
mostly by reading source definitions and equality cases rather than by
searching for counterexamples, which is also where most of the false starts
live.

**These results are now under external review, one conjecture at a time.**
The first round of feedback (29 Jul 2026) corrected conjecture 65 — a
variable had been misread, and the "trivial" proof was of a different,
easier statement. That correction is reflected above and in the log; more
feedback will land the same way. Until a conjecture has been confirmed
through that review, treat every "proved" claim here as this project's
reading of the formal statement, not as independently verified.

## Where to look

- **Full narrative, including every failed attempt:** [RESEARCH_LOG.md](RESEARCH_LOG.md)
- **File-by-file index of the code:** the "Files" section at the end of the log
- **Requirements:** Python 3, `numpy`, `networkx`, `mpmath`
