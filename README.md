# BMSSP — C++ Implementation of "Breaking the Sorting Barrier for Directed SSSP"

A C++ implementation and empirical benchmark of the **Bounded Multi-Source Shortest Path (BMSSP)** algorithm from:

> Ran Duan, Jiayi Mao, Xiao Mao, Xinkai Shu, Longhui Yin.
> *Breaking the Sorting Barrier for Directed Single-Source Shortest Paths.*
> arXiv:2504.17033, July 2025.

This paper gives the first deterministic algorithm that solves single-source shortest paths (SSSP) on directed, non-negatively weighted graphs in **O(m log^(2/3) n)** time in the comparison-addition model — asymptotically faster than Dijkstra's classic **O(m + n log n)** bound on sparse graphs, and the first result to show Dijkstra's algorithm is *not* optimal for the plain distance-only version of SSSP.

This repository contains a from-scratch reimplementation of the algorithm (`FindPivots`, `BaseCase`, `BMSSP`, and the partial-sorting structure `D`), plus a small driver that generates a random graph, runs both Dijkstra and BMSSP on it, and reports timing and correctness statistics.

> ⚠️ **This is an experimental / educational implementation**, not a faithful, asymptotically-optimal realization of every data structure in the paper. See [Known Simplifications & Deviations](#known-simplifications--deviations-from-the-paper) below before relying on it for correctness or for reproducing the paper's performance claims.

---

## Table of Contents

- [Background](#background)
- [Algorithm Overview](#algorithm-overview)
- [Repository Structure](#repository-structure)
- [Build Instructions](#build-instructions)
- [Usage](#usage)
- [Code ↔ Paper Mapping](#code--paper-mapping)
- [Known Simplifications & Deviations from the Paper](#known-simplifications--deviations-from-the-paper)
- [Interpreting the Benchmark Output](#interpreting-the-benchmark-output)
- [Limitations & Possible Improvements](#limitations--possible-improvements)
- [References](#references)
- [License](#license)

---

## Background

For decades, the fastest general SSSP algorithms have been built around one of two ideas:

- **Dijkstra's algorithm** — repeatedly extract the closest unvisited vertex via a priority queue and relax its outgoing edges. Because it fundamentally sorts vertices by distance, it costs at least Θ(n log n) in the comparison-addition model.
- **Bellman–Ford** — relax *all* edges for a fixed number of rounds. This avoids sorting entirely, but costs O(mk) time to correctly resolve paths with up to k edges.

The paper's key insight is that Dijkstra's bottleneck comes from having to fully order a potentially large "frontier" of active vertices. By recursively partitioning the vertex set into a hierarchy of bounded sub-problems, and only fully processing a *small* set of "pivot" vertices at each level (roughly `|frontier| / k` of them, where `k = log^(1/3) n`), the algorithm avoids ever sorting more than a shrinking sub-linear slice of the graph. Combining this pivoting trick with short bursts of Bellman-Ford-style relaxation yields an algorithm whose per-vertex overhead drops from `O(log n)` to `O(log n / log^Ω(1) n)`, i.e. **O(m log^(2/3) n)** total.

## Algorithm Overview

The algorithm is built from three pieces, called recursively over `O(log n / t)` levels (with `k = log^(1/3) n`, `t = log^(2/3) n`):

1. **`FindPivots(B, S)`** (Lemma 3.2 / Algorithm 1)
   Given a frontier `S` and bound `B`, relax edges from `S` for `k` rounds. If the explored region `W` stays small (`|W| ≤ k|S|`), collapse `S` down to a small set of "pivot" roots `P` (each pivot corresponds to a shortest-path subtree of size ≥ k), so `|P| ≤ |W| / k`. Otherwise fall back to `P = S`.

2. **`BaseCase(B, S)`** (Algorithm 2)
   When the recursion bottoms out (`l = 0`, `S` a single complete vertex), run a small bounded Dijkstra from that vertex to find its `k` closest neighbors within distance `B`.

3. **`BMSSP(l, B, S)`** (Algorithm 3)
   The main recursive driver. It calls `FindPivots` to shrink `S` to pivots `P`, inserts `P` into a specialized partial-sorting data structure `D` (Lemma 3.3), then repeatedly:
   - Pulls a small batch `(Bᵢ, Sᵢ)` of the currently-smallest keys out of `D`,
   - Recurses one level down: `BMSSP(l-1, Bᵢ, Sᵢ)`,
   - Relaxes edges out of the returned complete set, feeding newly-discovered vertices back into `D` (either via `Insert` or via a `Batch Prepend`, depending on whether their new distance falls before or after the current pull boundary),
   - Stops early if the accumulated output set `U` grows past `Θ(k·2^(lt))`, signaling a "partial execution" that reports a smaller boundary `B' < B`.

   Data structure `D` (Lemma 3.3) is the piece that lets pivots be inserted/batch-prepended and pulled out again without ever fully sorting them — its `Insert`, `Batch Prepend`, and `Pull` operations run in `O(log(N/M))`, amortized `O(log(L/M))` per element, and `O(|S'|)` respectively.

The top-level call is `BMSSP(⌈log n / t⌉, ∞, {s})`, which (because the output set can be at most all of `V`) always terminates as a "successful execution" and returns exact distances to every reachable vertex.

## Repository Structure

```
.
├── prj.txt                 # BMSSP implementation + Dijkstra benchmark driver (C++)
└── README.md                # This file
```

> The source currently lives in `prj.txt`. It is plain C++17 — rename it to e.g. `bmssp.cpp` (or update the build command below) if you'd like it to compile with standard tooling / IDE syntax highlighting.

Main components inside the source file:

| Symbol | Role |
|---|---|
| `struct Edge` | `(dest, weight)` directed edge |
| `struct Pivot` | Output of `FindPivot`: pivot set `P`, explored frontier `W` |
| `struct Result` | Output of `BaseCase` / `BMSSP`: boundary `B` and complete vertex set `U` |
| `class DataStructureD` | Simplified stand-in for the Lemma 3.3 block-linked-list structure, backed by two `std::priority_queue`s + lazy deletion |
| `FindPivot(...)` | Implements (a simplified version of) Algorithm 1 |
| `BaseCase(...)` | Implements Algorithm 2 (bounded mini-Dijkstra) |
| `BMSSP(...)` | Implements Algorithm 3 (the main recursive procedure) |
| `dijkstra(...)` | Textbook Dijkstra, used as the correctness/performance baseline |
| `main()` | Random graph generator + benchmark driver |

## Build Instructions

Requires a C++17-capable compiler (uses structured bindings, `#include "bits/stdc++.h"`).

```bash
# rename or copy the source to a .cpp file first
cp prj.txt bmssp.cpp

g++ -O2 -std=c++17 -o bmssp bmssp.cpp
```

## Usage

The program reads three integers from **stdin**: `n m seed`.

```bash
./bmssp
# then type, e.g.:
# 100000 300000 42
```

or pipe input directly:

```bash
echo "100000 300000 42" | ./bmssp
```

**Input parameters**

| Param | Meaning |
|---|---|
| `n` | Number of vertices |
| `m` | Target edge count *(currently unused — see [deviations](#known-simplifications--deviations-from-the-paper))* |
| `seed` | RNG seed for reproducible graph generation |

**Output**

```
Dijkstra time: <seconds>
BMSSP time   : <seconds>
Max diff     : <largest |d_dijkstra(v) - d_bmssp(v)| over commonly-reached v>
Common reachable: <number of vertices both algorithms found reachable>
```

- `Max diff` should be `0` if BMSSP computed exactly the same distances as Dijkstra for every vertex both algorithms reached — this is the correctness check.
- `Common reachable` should equal `n` (every vertex is reachable, since the generated graph is a spanning tree rooted at vertex 0) if both algorithms fully explored the graph.

## Code ↔ Paper Mapping

| Paper | Code |
|---|---|
| Lemma 3.2 / Algorithm 1 (`FindPivots`) | `FindPivot(int k, edges, B, S, dist, p_limit)` |
| Algorithm 2 (`BaseCase`) | `BaseCase(int k, edges, B, S, dist)` |
| Lemma 3.3 (partial-sorting structure `D`) | `class DataStructureD` (`insert`, `Batch_Prepend`, `pull`, `isEmpty`) |
| Algorithm 3 (`BMSSP`) | `Result BMSSP(int l, edges, B, S, dist)` |
| `k = ⌊log^(1/3) n⌋`, `t = ⌊log^(2/3) n⌋` | `int k = max(1,(int)pow(log2(n),1.0/3.0))`, `int t = max(1,(int)pow(log2(n),2.0/3.0))` inside `BMSSP` |
| `M = 2^((l-1)t)` (block size for `D`) | `int M = (int) pow(2, max(0,(l-1)*t));` |
| Top-level call `BMSSP(⌈log n/t⌉, ∞, {s})` | `BMSSP(l, adj, LLONG_MAX, {src}, d_bm);` in `main()` |

## Known Simplifications & Deviations from the Paper

This implementation prioritizes being a runnable, testable approximation of the algorithm's *logical structure* over reproducing every asymptotic guarantee from the paper. If you're comparing this code against the paper line-by-line, note the following:

1. **`DataStructureD` is not the Lemma 3.3 block-linked-list structure.**
   The paper's structure achieves amortized `O(max{1, log(N/M)})` inserts and `O(L·max{1, log(L/M)})` batch-prepends via a two-sequence, self-balancing block layout. This code instead uses two ordinary `std::priority_queue`s (one for `Insert`, one for `Batch Prepend`) with lazy deletion via a `best[]` hash map. This is simple and correct in spirit (still respects "smallest value per key", still separates freshly-batch-prepended keys from inserted ones during `pull`), but every operation costs `O(log n)`, not the paper's finer amortized bounds — so the code will **not** exhibit the paper's `O(m log^(2/3) n)` scaling in practice, even though the recursive algorithmic structure is preserved. The `Batch_Prepend` method's chunking loop (splitting into blocks of size `M/2`) has no effect on behavior since everything still lands in one shared heap — it's a structural echo of the paper's block-splitting, not a functioning equivalent.

2. **`FindPivot` does not implement the exact tree-based pivot rule.**
   Algorithm 1 in the paper relaxes edges for `k` rounds while updating `db[·]`, then builds the shortest-path forest restricted to `W` and selects roots whose subtree has ≥ `k` vertices as pivots. In this code:
   - The `k`-round frontier expansion (building `W`) is present, but it **only checks** `dist[u] + weight < B` — it does not write the relaxed value back into `dist[]`. So, unlike the paper, exploring `W` here does not make any vertex "complete."
   - The forest-construction / subtree-size pivot selection (the commented-out block near the end of `FindPivot`) is disabled. Instead, pivots `P` are chosen heuristically: the `p_limit` closest candidate vertices in `S` (by current `dist`), or an arbitrary prefix of `S` if none qualify.
   - `p_limit` itself is capped at `2^min(10, t)` rather than the paper's uncapped `2^t`, to keep pivot-set sizes practical for large `n`.

   This makes `FindPivot` a heuristic approximation rather than a certified `O(min{k²|S|, k|Uₑ|})`-time routine with the `|P| ≤ |W|/k` guarantee from Lemma 3.2.

3. **`B₀` computation in `BMSSP` has a non-paper guard clause.**
   Lemma 3.1 / Algorithm 3 sets `B'₀ ← min_{x∈P} db[x]` unconditionally (or `B` if `P` is empty). The code instead does:
   ```cpp
   if (!pivot.P.empty() and *max_element(pivot.P.begin(), pivot.P.end()) != 0) {
       for (int x : pivot.P) B0 = min(B0, dist[x]);
   }
   ```
   The extra `max_element(...) != 0` condition (skipping the min-update whenever pivot vertex `0` — i.e. the source — is the largest-indexed pivot) is not part of the paper's algorithm and looks like an ad-hoc guard rather than an intentional optimization. It can leave `B0` at its stale value (`B`) in some cases where the paper would have tightened it.

4. **The benchmark graph is a random tree, not a general sparse graph.**
   In `main()`, only the "Backbone" loop is active — each vertex `i` gets a single random parent among `0..i-1`, producing `n - 1` edges total (a random recursive tree). The "Remaining edges" loop that would add extra random directed edges up to the requested `m` is present in the source but commented out, so **the `m` input parameter is currently unused**, and the generated graph has no cycles, no parallel paths, and no vertices with in-degree/out-degree > O(1). This is a much easier instance than the general constant-degree directed graphs the paper's complexity analysis targets, so timing comparisons on this generator will not reflect worst-case behavior. Uncomment / adapt that loop (and feed `m` into it) if you want denser test graphs.

5. **No constant-degree graph transformation.**
   The paper assumes (Section 2, "Constant-Degree Graph") that the input graph has been transformed so every vertex has O(1) in/out-degree, using the cycle-substitution construction attributed to Frederickson. This code runs directly on the input graph as given (which happens to already be constant-degree here, since it's a random tree with a single parent edge per vertex) — the transformation itself is not implemented.

6. **`k`/`t` in `BMSSP` vs. `k`/`t`/`l` in `main()` are computed differently.**
   Inside `BMSSP`, `k` and `t` are derived from `log2(n)`. In `main()`, the top-level recursion depth `l` is derived independently from `log(n)` (natural log) and a locally-recomputed `t`. These two computations are not kept in sync via shared constants — they happen to be asymptotically consistent, but a refactor extracting `k`/`t` into shared, single-source parameters would make the relationship clearer and less error-prone.

7. **`Assumption 2.1` (unique path lengths / total tie-breaking order) is not enforced.**
   The paper handles ties between equal-length paths via a lexicographic tuple order over `(length, hop-count, vertex sequence)` so that predecessor pointers stay well-defined. This code uses plain `long long` distance comparisons with `<=` for relaxation (see `BaseCase`, `BMSSP`, `dijkstra`), so ties are broken arbitrarily by insertion/processing order rather than by the paper's canonical rule. This has no effect on final shortest-distance values (which is all this benchmark checks), but would matter if you needed a deterministic/reproducible shortest-path *tree*.

None of the above invalidate the algorithm's recursive logic (`FindPivots → pivots → recursive BMSSP calls → relax → batch-prepend`) — the control flow genuinely mirrors Algorithm 3. But the *complexity* guarantees and a couple of correctness-adjacent details (points 2 and 3 above, in particular) diverge from the paper, which is worth knowing before trusting `Max diff == 0` as a proof of general correctness, or before citing this code's runtime as evidence of `O(m log^(2/3) n)` behavior.

## Interpreting the Benchmark Output

- **`Max diff == 0` and `Common reachable == n`**: BMSSP agreed with Dijkstra on every vertex's distance, on this particular random-tree instance.
- **`Max diff > 0`**: BMSSP diverged from the true shortest-path distances on at least one vertex for this run — likely a consequence of one of the simplifications above (most plausibly #2, since pivots/exploration don't always behave exactly as the paper specifies) rather than a fundamental flaw in the recursive skeleton. Worth investigating with a smaller `n` and the commented-out debug `cout` statements re-enabled.
- **`Common reachable < n`**: some vertices were not marked reachable by one or both algorithms — on the tree generator this should not happen for Dijkstra; if it happens for BMSSP specifically, it indicates the recursion terminated (successful execution) without actually completing all vertices, which the paper's Lemma 3.9 guarantees shouldn't occur for a top-level call with `B = ∞`.
- **Timing**: because `DataStructureD` uses generic `O(log n)` heaps rather than the paper's amortized structure, and because the test graph is a simple random tree (Dijkstra is already near-optimal on trees), do not expect to see BMSSP outperform Dijkstra in wall-clock time on this benchmark. The value of this code is in exercising the algorithm's logical structure and checking correctness, not in demonstrating the paper's asymptotic speedup.

## Limitations & Possible Improvements

- [ ] Implement the true Lemma 3.3 block-linked-list structure (or another structure with matching amortized bounds) in place of the dual-`priority_queue` `DataStructureD`.
- [ ] Restore the forest-based pivot selection in `FindPivot` (including writing back relaxed distances during the k-round expansion) to match Algorithm 1 exactly, and remove the `p_limit` cap.
- [ ] Investigate and fix/remove the `max_element(...) != 0` guard in `BMSSP`'s `B0` computation.
- [ ] Re-enable / parameterize the "remaining edges" generator in `main()` so `m` is actually used and denser, cyclic graphs can be benchmarked.
- [ ] Add the constant-degree graph transformation described in Section 2 of the paper for general input graphs.
- [ ] Add unit tests on small hand-checkable graphs (not just statistical comparison against Dijkstra on random trees).
- [ ] Extract `k`, `t` computation into a single shared helper used by both `main()` and `BMSSP()`.

## References

- Ran Duan, Jiayi Mao, Xiao Mao, Xinkai Shu, Longhui Yin. *Breaking the Sorting Barrier for Directed Single-Source Shortest Paths.* arXiv:2504.17033v2 [cs.DS], July 2025.
- E. W. Dijkstra. *A note on two problems in connexion with graphs.* Numerische Mathematik, 1959.
- R. Bellman. *On a routing problem.* Quarterly of Applied Mathematics, 1958.

## License

No license file is currently included in this repository. Add one (e.g. MIT, Apache-2.0) if you intend for others to reuse this code.
