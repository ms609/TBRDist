# TBRDist — Agent Memory

## Package overview
R package for phylogenetic tree rearrangement distances. Wraps two C++ libraries
(uspr and rspr) as git submodules in `src/uspr/` and `src/rspr/`.
All C++ is C++17; both submodules compile as single translation units
(`src/uspr.cpp`, `src/rspr.cpp`) to avoid multiple-definition errors from
header-only globals. Exports `TBRDist()`, `USPRDist()`, `ReplugDist()`
(unrooted, via uspr) and `RSPRDist()` (rooted, via rspr).

## Version
Currently **1.1.0** (ready for CRAN). All 147 tests pass.
`R CMD check` gives 0 errors, 0 warnings, 1 NOTE (R 4.5 bug in
`read_symbols_from_dll` — unrelated to TBRDist, documented in
`cran-comments.md`).

## Dependencies
- `TreeTools (>= 2.1.0.9003)` — required because `tree_number.h` is defined
  in `TreeTools/inst/include/TreeTools/tree_number.h`. This header was written
  for TBRDist but lives in TreeTools as its natural home (it implements the C++
  backend of TreeTools' tree-numbering functionality). TreeTools 2.1.0.9003 is
  the first version containing it; pin will be relaxed once TreeTools releases
  to CRAN.
- `LinkingTo: BH, Rcpp, TreeTools`

## Architecture: numbered A* search (the fast path)
After the initial leaf-reduction, if the reduced tree has ≤ 51 leaves the
search uses `uspr_distance_numbered()` in `src/uspr/uspr.h`.

### Key data structures
- **Tree identity**: 256-bit `tree_num_t` (4 × uint64_t, Tromp encoding),
  defined in `TreeTools/inst/include/TreeTools/tree_number.h`.
- **Visited set**: `unordered_set<tree_num_t, tree_num_hash>` — O(1) lookup.
- **Priority queue**: `multiset<tree_distance_num>` — sorted by
  (total_distance, estimate, estimator).
- **Neighbor generation**: `get_neighbors_numbered()` in
  `src/uspr/uspr_neighbors_numbered.h` — applies each SPR, encodes result as
  `tree_num_t`, reverts. Never serializes to Newick.

### Decode path (each queue pop)
```
tree_num_t → tree_number_to_parent() → direct utree construction
           → uforest(utree&&) move ctor → normalize_order()
```
`tree_number_to_utree()` builds the utree directly from the parent vector (no
Newick). The R root node (2n−1) is dissolved — its two children are connected
directly, producing a proper unrooted degree-3 tree with exactly n−2 internal
nodes. Defined in `src/uspr/tree_numbering.h`.

### Estimator cascade (progressive A*)
Each unique tree number is popped up to 4× with escalating estimators:
**BFS → TBR_APPROX → TBR → REPLUG**

This cascade is **intentional design**: TBR is much cheaper to compute than
REPLUG, so calling TBR first provides a tighter bound cheaply before
committing to the expensive REPLUG computation. Do NOT remove the TBR level.

Each estimator function (`tbr_high_lower_bound`, `tbr_distance`,
`replug_distance` in `src/uspr/tbr.h`) makes its own deep copy and calls
`root()` + `distances_from_leaf_decorator` internally.

### Phase 2/3 lookup tables
- **Phase 2**: After leaf_reduction, if 4–9 leaves remain, look up exact SPR
  distance from precomputed tables (`src/spr/` + `src/spr_lookup.cpp`).
- **Phase 3 (BFS-level only)**: On first pop of each tree (estimator == BFS),
  run a mini leaf-reduction and attempt table lookup on the result (4–9
  leaves). Fires on ~50% of BFS-level pops.

### `distances_from_leaf_decorator` — NOT needed in numbered path
`distances_from_leaf_decorator(T, leaf)` sets `distance` on every node via BFS
from `leaf`. It is used by `cut_edge()` and `contract()` in the TBR forest
algorithms. It is **not** needed by:
- `normalize_order()` (uses leaf labels recursively, not distances)
- `utree_to_tree_number()` (postorder DFS, no distances)
- `utree::uspr()` (neighbor-list manipulation, no distances)

The outer A* loop therefore does NOT call it; the Phase 3 copies do call it
(feeding into `leaf_reduction_hlpr` which uses `cut_edge`).

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
- Neighbor lists are **bidirectional** (unrooted). `get_parent()` =
  `neighbors.front()`.
- After `normalize_order()`, each node's neighbor list is sorted: parent
  first, then children by smallest-descendant leaf label.

## Benchmark results (same-session A/B, 15 pairs, 10–12 leaves, Windows)
| Version | Mean | Notes |
|---------|------|-------|
| Newick roundtrip in A* loop | 11.68 s | baseline |
| Direct utree construction | 11.18 s | −4.3% |
| + removed redundant `distances_from_leaf_decorator` | 11.30 s | ≈0% (noise) |
| + `priority_queue` (binary min-heap) replacing `multiset` | 8.94 s | **−21%** |

**Key finding**: The `multiset` red-black tree caused significant cache-miss
overhead. Replacing with `priority_queue` (contiguous `vector` storage) gave
the largest single speedup. The O(n) per-pop items are negligible; queue
operations were the hidden cost.

## Known remaining optimisation candidates
1. **`multiset` priority queue** — O(log Q) insert/erase with poor cache
   behaviour (scattered red-black tree nodes). Replacing with `priority_queue`
   (binary min-heap, contiguous storage) keeps O(log Q) complexity but
   improves cache performance. Worth trying; the multiset is not known to be
   a bottleneck but hasn't been measured either.
2. **Phase 3 per-pop overhead** — every BFS pop deep-copies T and T2 and runs
   `leaf_reduction_hlpr`. Could be gated on n > some threshold.
3. **`normalize_order()` constant factor** — uses `map<int, unode*>` per
   internal node. Each internal node has exactly 2 children; the map could be
   replaced with a two-element compare-and-swap for a small constant-factor win.
4. **Memory allocation patterns** — `list<unode*>` for neighbor lists causes
   many small heap allocations. At scale, fixed-size arrays or an arena
   allocator would improve cache behaviour.
5. **Extending lookup tables to 10–11 leaves** — would eliminate A* for pairs
   reducing to ≤11 leaves. Tables too large to ship in the package per user.
