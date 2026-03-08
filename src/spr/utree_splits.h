#pragma once

// Bridge between uspr's utree representation and the splits-based lookup tables.
// Converts a utree to a SplitSet for exact SPR distance lookup on small trees.

#include <cstdint>
#include <array>
#include "lookup.h"
#include "../uspr/unode.h"
#include "../uspr/utree.h"

// Defined in spr_lookup.cpp
int lookup6(const SplitSet6& sp1, const SplitSet6& sp2);
int lookup7(const SplitSet7& sp1, const SplitSet7& sp2);
int lookup8(const SplitSet8& sp1, const SplitSet8& sp2);
int lookup9(const SplitSet9& sp1, const SplitSet9& sp2);

namespace spr_lookup {

// DFS collecting descendant-leaf bitmasks.
// Each internal node (except the root internal node) contributes one split.
template <typename SplitT>
SplitT collect_splits_dfs(const unode* node, const unode* parent,
                          SplitT* splits, int& split_idx,
                          bool is_root_internal) {
  if (node->get_label() >= 0) {
    return static_cast<SplitT>(1) << node->get_label();
  }

  SplitT mask = 0;
  for (unode* child : node->const_neighbors()) {
    if (child != parent) {
      mask |= collect_splits_dfs<SplitT>(child, node, splits, split_idx, false);
    }
  }

  if (!is_root_internal) {
    splits[split_idx++] = mask;
  }

  return mask;
}

// Extract the n-3 non-trivial splits from a utree into a SplitSet.
// Leaves must be labelled 0..n-1.
template <typename SplitSetT, typename SplitT>
bool utree_to_splits(utree& T, SplitSetT& out) {
  int sl = T.get_smallest_leaf();
  unode* root_leaf = T.get_leaf(sl);
  if (root_leaf == nullptr) return false;

  const auto& nbrs = root_leaf->const_neighbors();
  if (nbrs.empty()) return false;
  unode* root_internal = nbrs.front();

  SplitT buf[8]; // max 6 splits (9 leaves)
  int idx = 0;

  collect_splits_dfs<SplitT>(root_internal, root_leaf, buf, idx, true);

  constexpr int expected = std::tuple_size<SplitSetT>::value;
  if (idx != expected) return false;

  for (int i = 0; i < expected; ++i) {
    out[i] = buf[i];
  }
  return true;
}

inline int count_leaves(utree& T) {
  int n = 0;
  for (unode* u : T.get_leaves()) {
    if (u != nullptr) ++n;
  }
  return n;
}

// Attempt exact SPR lookup for trees with 6-9 leaves.
// Returns the exact distance, or -1 if trees are too large/small.
inline int lookup_utrees(utree& T1, utree& T2) {
  int n = count_leaves(T1);

  switch (n) {
  case 6: {
    SplitSet6 sp1{}, sp2{};
    if (!utree_to_splits<SplitSet6, Split6>(T1, sp1)) return -1;
    if (!utree_to_splits<SplitSet6, Split6>(T2, sp2)) return -1;
    // lookup6 handles smaller_split internally
    return lookup6(sp1, sp2);
  }
  case 7: {
    SplitSet7 sp1{}, sp2{};
    if (!utree_to_splits<SplitSet7, Split7>(T1, sp1)) return -1;
    if (!utree_to_splits<SplitSet7, Split7>(T2, sp2)) return -1;
    for (auto& s : sp1) s = smaller_split7(s);
    return lookup7(sp1, sp2);
  }
  case 8: {
    SplitSet8 sp1{}, sp2{};
    if (!utree_to_splits<SplitSet8, Split8>(T1, sp1)) return -1;
    if (!utree_to_splits<SplitSet8, Split8>(T2, sp2)) return -1;
    for (auto& s : sp1) s = smaller_split8(s);
    return lookup8(sp1, sp2);
  }
  case 9: {
    SplitSet9 sp1{}, sp2{};
    if (!utree_to_splits<SplitSet9, Split9>(T1, sp1)) return -1;
    if (!utree_to_splits<SplitSet9, Split9>(T2, sp2)) return -1;
    for (auto& s : sp1) s = smaller_split9(s);
    return lookup9(sp1, sp2);
  }
  default:
    return -1;
  }
}

} // namespace spr_lookup
