# TBRDist — Agent Memory

## Standing instructions

**Always obtain human approval before running `git push` or any command that
writes to a remote.** Stage and commit freely, but show the diff/commit summary
and wait for an explicit go-ahead before pushing.

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
| `src/uspr/utree.h` | Default ctor, move ctor, arena-based allocation |
| `src/uspr/unode.h` | `neighbor_list`, `contracted_list` inline arrays |
| `src/uspr/unode_arena.h` | Block arena allocator for unode objects |
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
- **Memory management**: All unodes are allocated from a `unode_arena` owned by
  each `utree`. The arena allocates in blocks of 32 unodes (placement new into
  raw memory). Deallocation is bulk — the arena frees entire blocks at once.
  `unode` is trivially destructible so no per-object teardown is needed.
  `add_internal_node()`, `add_leaf()`, and the copy constructor all use the arena.

## Allocation optimizations applied to tbr.h branch-and-bound
The recursive `tbr_distance_hlpr` in `src/uspr/tbr.h` deep-copies F1, F2, and
bookkeeping structures at each branch (up to 5 branches per call, exponential
depth). The following changes reduced per-copy overhead:

- **`neighbor_list`**: `list<unode*>` → fixed inline array (capacity 4) in
  `unode.h`. Zero heap allocation per node.
- **`contracted_list`**: `vector<unode*>` → fixed inline array (capacity 3).
- **`unode_arena`**: Block arena allocator (32 unodes/block) replaces individual
  `new`/`delete` for every node in every tree copy.
- **`singletons`**: `list<int>` → `vector<int>`.
- **`nodemapping`**: `map<int,int>` → `unordered_map<int,int>` for forward/
  backward mappings (point lookups only, never iterated).
- **`find_leaves()`**: returns `vector<int>` instead of `list<int>`.
- **`find_pendants()`**: returns `vector<pair<int,int>>`.
- **`get_path()`**: uses `vector<unode*>` internally.
- **`sibling_pairs`**: kept as `map<int,int>` (iteration order matters for
  deterministic MAF output that tests depend on).

## Benchmark methodology
**Fixed pairs** (reproducible): 60 random 10-tip pairs generated via
`as.phylo(floor(runif(60)*NUnrooted(10L)), nTip=10L)` with `set.seed(3141)`.
All pairs have SPR distances 3–5 (table: 5×d=3, 40×d=4, 15×d=5). Using
nearby-index (small-offset) pairs is NOT appropriate: distance 1–2 pairs are
leaf-reduced to ≤9 leaves and go straight to the lookup table, never
exercising the A* loop or `normalize_order`.

## Benchmark results
### Session 1 — same-session A/B, 15 pairs, 10–12 leaves, Windows
| Version | Time | Notes |
|---------|------|-------|
| Newick roundtrip in A* loop | 11.68 s | baseline |
| Direct utree construction | 11.18 s | −4.3% |
| + removed redundant `distances_from_leaf_decorator` | 11.30 s | ≈0% (noise) |
| + `priority_queue` (binary min-heap) replacing `multiset` | 8.94 s | **−21%** |

### Session 2 — fixed 60-pair benchmark, 10-tip random, Windows
| Version | Time | Notes |
|---------|------|-------|
| `map`-based `normalize_order_hlpr` | 19.14 s | baseline |
| inline 3-array + insertion sort | 17.88–19.31 s | ≈0% (noise) |

### Session 3 — allocation optimizations, 60-pair benchmark, Windows
| Version | Time | Δ from prev |
|---------|------|-------------|
| Inline neighbor/contracted arrays | 9.72 s | −49% from 19.14 s |
| + container replacements (vectors, unordered_maps) | 9.94 s | ≈0% (noise) |
| + unode arena allocator (block size 32) | 7.20 s | **−26%** |

**Cumulative**: 19.14 s → 7.20 s = **2.66× faster** (62% reduction).

### VTune profile (pre-arena, at 9.72 s)
- 63% of runtime in `malloc_base` (38%) + `free_base` (25%)
- Top C++ hotspots: `utree::utree` copy ctor (0.57 s),
  `vector::_M_realloc_append` for contracted_neighbors (0.33 s),
  `list<int>::_M_clear` for singletons (0.14 s)
- TBR/replug branch-and-bound dominates: exponential copies amplify all
  allocation costs

### VTune profile (post-arena, at 7.20 s / 9.60 s CPU)
- 54% of runtime in `malloc_base` (31%) + `free_base` (23%) — down from 63%
- Top hotspots in TBRDist.dll:
  - `utree::utree` copy ctor: 0.65 s (6.7%) — mostly container allocs now,
    not per-node new/delete
  - `func@0x1c5f79f70`: 0.32 s (3.3%) — unresolved, likely a TBR helper
  - `vector<unode*>::_M_realloc_append`: 0.11 s — `internal_nodes`/`leaves`
    vectors inside each utree copy
  - `tbr_distance_hlpr`: 0.11 s — the branch-and-bound logic itself
  - `tbr_approx_hlpr`: 0.10 s
  - `unode::contract`: 0.09 s
  - `uforest::cut_edge`: 0.06 s
  - `nodemapping::~nodemapping`: 0.03 s — unordered_map destruction
  - `unordered_map::_M_insert_unique_node`: 0.05 s
  - `map::insert` (sibling_pairs): 0.03 s
- Remaining malloc/free is spread across many small sources (vector/map/
  unordered_map containers inside each utree copy) rather than concentrated
  in one dominant source.

**Key findings**:
- The `multiset` → `priority_queue` change was the dominant early win (−21%).
- Inline neighbor arrays produced a large win (−49%) because the copy
  constructor allocates per-node neighbor storage — eliminating those
  allocations directly reduced the dominant malloc/free cost.
- Container replacements (vectors, unordered_maps) were individually
  unmeasurable but are kept for code quality.
- The arena allocator produced a further −26% by batching unode allocations
  into 32-node blocks, reducing malloc/free calls by ~32×.
- Arena block size does not matter: 32, 64, 128 all benchmark within noise
  (7.20–7.55 s).

## Known remaining optimisation candidates
1. **Phase 3 per-pop overhead** — every BFS pop deep-copies T and T2 and runs
   `leaf_reduction_hlpr`. Could be gated on n > some threshold to avoid the
   copy overhead when the reduced tree cannot possibly hit the 4–9-leaf window.
2. **TBR/replug branch-and-bound** — the measured dominant cost. Speedups here
   require either better pruning (algorithmic, hard) or profiling to find any
   unexpected constant-factor bottleneck within `tbr.h`.
3. **Extending lookup tables to 10–11 leaves** — would eliminate A* for pairs
   reducing to ≤11 leaves. Tables too large to ship in the package per user.
4. **Revert-based SPR in branch-and-bound** — modify trees in place and undo
   after recursing, instead of deep-copying. Saves ~1 copy per branch point
   (first branch only), not all copies. High implementation risk for moderate gain.
5. **Identify unresolved `func@0x1c5f79f70`** — 3.3% of CPU time. Building
   with full debug symbols may resolve this to a specific function.
