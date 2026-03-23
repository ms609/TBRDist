# TBRDist 2.0.0 (2026-03-23)

- Add `RSPRDist()` for rooted Subtree Prune-and-Regraft (rSPR) distances,
  wrapping the exact FPT algorithm (with cluster decomposition) of Whidden,
  Beiko & Zeh (2013). Supports exact distances, a linear-time 3-approximation,
  and maximum agreement forest output.
- Improve performance of `USPRDist()` by replacing string-based tree
  representation in the A\* search with 256-bit integer tree numbers
  (Tromp encoding), giving O(1) equality and hash lookup.
- Add exact SPR distance lookup tables for trees with 4–9 leaves, eliminating
  A\* search for small trees and subtrees encountered during reduction.
- Improve performance of `USPRDist()` A\* search by replacing the
  `multiset`-based priority queue (pointer-chasing red-black tree) with a
  binary min-heap (`priority_queue` over a contiguous `vector`), giving better
  cache behaviour for queue operations.
- Reduce heap allocations in `normalize_order()` by replacing per-node
  `map<int,unode*>` with a 3-element inline array sort.
- Reduce heap allocations in TBR/replug branch-and-bound search by replacing
  `list<unode*>` neighbour storage with fixed-capacity inline arrays, replacing
  `list`/`map` book-keeping containers with `vector`/`unordered_map`, and
  allocating tree nodes from a block arena instead of individual `new`/`delete`.

# TBRDist 1.0.3 (2025-11-28)

- Match updated CRAN policy.

# TBRDist 1.0.2

- Import RdMacros package 'Rdpack'.

# TBRDist 1.0.1

- Address memory mismanagement in `USPRDist()`.

# TBRDist 1.0.0

- Initial implementation of distances on unrooted trees.
