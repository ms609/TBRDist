# TBRDist — Agent Memory

## Package overview
R package for phylogenetic tree rearrangement distances. Wraps two C++ libraries
(uspr and rspr) as git submodules in `src/uspr/` and `src/rspr/`.
All C++ is C++17; both submodules are compiled as single translation units
(`src/uspr.cpp`, `src/rspr.cpp`) to avoid multiple-definition errors from
header-only globals. The package exports `TBRDist()`, `USPRDist()`, `ReplugDist()`
(unrooted, via uspr) and `RSPRDist()` (rooted, via rspr).

## Version
Currently **1.1.0** (ready for CRAN). All 147 tests pass.
`R CMD check` gives 0 errors, 0 warnings, 1 NOTE (R 4.5 bug in `read_symbols_from_dll`
— unrelated to TBRDist, documented in `cran-comments.md`).

## Architecture: numbered A* search (the fast path)
After the initial leaf-reduction, if the reduced tree has ≤ 51 leaves the
search uses `uspr_distance_numbered()` in `src/uspr/uspr.h`.
Key data structures:
- **Tree identity**: 256-bit `tree_num_t` (4 × uint64_t, Tromp encoding).
  Defined in `TreeTools/inst/include/TreeTools/tree_number.h`.
- **Visited set**: `unordered_set<tree_num_t, tree_num_hash>` — O(1) lookup.
- **Priority queue**: `multiset<tree_distance_num>` — sorted by cost+estimate.
- **Neighbor generation**: `get_neighbors_numbered()` in
  `src/uspr/uspr_neighbors_numbered.h` — applies each SPR, encodes result as
  tree_num_t, reverts. Never serializes to Newick.

### Decode path (each queue pop)
```
tree_num_t → tree_number_to_parent() → direct utree construction
           → uforest(utree&&) move ctor → [optional] normalize_order()
```
`tree_number_to_utree()` builds the utree directly from the parent vector
(no Newick). The R root node (label 2n−1) is dissolved — its two children
are connected directly. Defined in `src/uspr/tree_numbering.h`.

### Estimator cascade
Each tree number is popped up to 4× with escalating estimators
(BFS → TBR_APPROX → TBR → REPLUG). Each estimator function (in `src/uspr/tbr.h`)
makes its own deep copy and calls `root()` + `distances_from_leaf_decorator`
internally — so the outer A* loop does NOT need to maintain distance state.

### Phase 2/3 lookup tables
- **Phase 2**: After leaf_reduction, if 4–9 leaves remain, look up exact SPR
  distance from precomputed tables (`src/spr/` + `src/spr_lookup.cpp`).
- **Phase 3 (BFS-level only)**: On first pop of each tree (estimator == BFS),
  run a mini leaf-reduction and attempt table lookup on the result (4–9 leaves).
  This fires on ~50% of BFS-level pops.

## Key files modified / created
| File | Description |
|------|-------------|
| `src/uspr/uspr.h` | Main A* dispatch + numbered search loop |
| `src/uspr/uspr_neighbors_numbered.h` | SPR neighbor generation → tree numbers |
| `src/uspr/tree_numbering.h` | utree ↔ tree_num_t bridge |
| `src/uspr/utree.h` | Added default ctor + move ctor |
| `src/uspr/uforest.h` | Added `uforest(utree&&)` move ctor |
| `src/rspr.cpp` | Rcpp wrapper for rspr rooted-distance library |
| `R/rspr.R` | `RSPRDist()` user-facing function |
| `tests/testthat/test-rspr.R` | 33 regression tests for RSPRDist |
| `TreeTools/inst/include/TreeTools/tree_number.h` | 256-bit tree_num_t + encoding |

## utree / uforest internals
- Leaves labeled 0..n-1; internal nodes labeled -2, -3, ...
- `internal_nodes[-(label)-2]` accesses an internal node by label.
- Neighbor lists are **bidirectional** (unrooted). `get_parent()` = `neighbors.front()`.
- After `normalize_order()`, each node's neighbor list is sorted: parent first,
  then children by smallest-descendant leaf label.
- `distances_from_leaf_decorator(T, leaf)`: sets `distance` on every node via
  BFS from `leaf`. Used by `cut_edge()` (TBR forest operations) and
  `contract()` (edge contraction in forests). **NOT** needed by `normalize_order()`,
  `utree_to_tree_number()`, or `utree::uspr()`.

## Benchmark results (same-session A/B, 15 pairs, 10-12 leaves, Windows)
| Version | Mean | Notes |
|---------|------|-------|
| Newick roundtrip in A* loop | 11.68s | baseline |
| Direct utree construction (no Newick) | 11.18s | −4.3% |
| + no `distances_from_leaf_decorator` | 11.30s | ≈0% (noise) |

**Key finding**: At n=10–12, the O(n) decorator calls are negligible. The
dominant cost is the TBR/replug estimator computation (O(2^k × n)), which the
A* calls 3–4× per unique tree explored. No O(n) cleanup optimisation will
make a meaningful dent here.

## Known remaining bottlenecks (in priority order)
1. **Estimator cascade itself** — `tbr_distance` / `replug_distance` do
   branch-and-bound search O(n × 4^k). This is the true bottleneck. Reducing
   the number of unique trees that need full estimation (better pruning,
   better initial estimate) would help most.
2. **Skipping the TBR level** — cascade is BFS → TBR_APPROX → TBR → REPLUG
   (4 pops/tree). Removing the TBR level (→ 3 pops) saves 25% of estimator
   calls at the cost of slightly worse ordering. Worth benchmarking.
3. **Extending lookup tables to 10–11 leaves** — currently 4–9 leaves only.
   Would eliminate A* entirely for pairs that reduce to ≤11 leaves. User noted
   tables would be too large for the package at 10 leaves.
4. **`multiset` priority queue** — O(log Q) insert/erase. A bucket queue
   indexed by integer distance (values ≤ 20 in practice) would give O(1).
   Worth implementing if queue size becomes large.
5. **Phase 3 per-pop overhead** — every BFS pop deep-copies T and T2 and runs
   `leaf_reduction_hlpr`. Could be gated on n > some threshold where reduction
   is likely to succeed.
6. **`normalize_order()` per SPR** — O(n) with map overhead. Each internal
   node has exactly 2 children; could replace `map<int,unode*>` with a
   compare-swap-2 for meaningful constant-factor improvement.
7. **Memory allocation patterns** — `list<unode*>` for neighbor lists causes
   many small heap allocations. At scale, an arena allocator or fixed-size
   neighbour arrays would improve cache behaviour significantly.
