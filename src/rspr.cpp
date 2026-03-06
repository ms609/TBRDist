#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#include <climits>
#include <map>
#include <sstream>
#include <string>

// rspr is a header-only library: all function bodies and global variables are
// defined in its headers.  We include everything from this single translation
// unit to avoid multiple-definition errors from the header-defined globals.
// The same pattern is used for uspr in src/uspr.cpp.
#include "rspr/rspr.h"

// ── Helpers ───────────────────────────────────────────────────────────────────

// RAII streambuf swap: redirects std::cout to a sink on construction and
// restores it on destruction.  Prevents rspr's internal diagnostic output
// (452+ cout calls) from reaching the R console.
struct CoutSuppressor {
  std::ostringstream sink;
  std::streambuf *orig;
  CoutSuppressor() : orig(std::cout.rdbuf(sink.rdbuf())) {}
  ~CoutSuppressor() { std::cout.rdbuf(orig); }
};

// Reset rspr's process-level optimisation globals to the same defaults as
// rspr.cpp's DEFAULT_OPTIMIZATIONS=true + DEFAULT_ALGORITHM=true blocks.
// These persist between R-level calls within a session, so must be reset
// before each use.
static void rspr_set_defaults() {
  // Algorithm
  BB                          = true;
  PREFER_RHO                  = true;
  // Branching optimisations
  CUT_ALL_B                   = true;
  CUT_ONE_B                   = true;
  REVERSE_CUT_ONE_B           = true;
  REVERSE_CUT_ONE_B_3         = true;
  CUT_TWO_B                   = true;
  CUT_AC_SEPARATE_COMPONENTS  = true;
  EDGE_PROTECTION             = true;
  EDGE_PROTECTION_TWO_B       = true;
  // Search order
  NEAR_PREORDER_SIBLING_PAIRS = true;
  PREORDER_SIBLING_PAIRS      = true;
  DEEPEST_PROTECTED_ORDER     = true;
  DEEPEST_ORDER               = true;
  PREFER_NONBRANCHING         = true;
  // Reductions
  LEAF_REDUCTION              = true;
  LEAF_REDUCTION2             = true;
  // Approximation optimisations
  APPROX_CUT_ONE_B            = true;
  APPROX_CUT_TWO_B            = true;
  APPROX_REVERSE_CUT_ONE_B    = true;
  // Clustering
  CLUSTER_TUNE                = 30;
  // Limits
  MAX_SPR                     = 1000;
  CLUSTER_MAX_SPR             = MAX_SPR;
  // Disabled by default
  MULTIFURCATING              = false;
  ALL_MAFS                    = false;
  VERBOSE                     = false;
}

// ── Exported function ─────────────────────────────────────────────────────────

//' Rooted SPR distance (C++ core)
//'
//' Low-level Rcpp interface to the rspr library.  Prefer the R wrapper
//' \code{RSPRDist()} for normal use.
//'
//' @param tree1,tree2 Character vectors of Newick-format rooted binary trees.
//'   Must be the same length.
//' @param approx Logical scalar: compute the linear-time 3-approximation?
//' @param exact  Logical scalar: compute the exact branch-and-bound distance?
//' @return A named list with elements \code{exact} (integer), \code{approx}
//'   (integer), \code{maf_1} and \code{maf_2} (character).  Elements not
//'   requested are filled with \code{NA} or empty strings.
//' @keywords internal
// [[Rcpp::export]]
List rspr_dist(const StringVector tree1,
               const StringVector tree2,
               const LogicalVector approx,
               const LogicalVector exact) {

  if (tree1.size() != tree2.size()) {
    throw std::length_error("Number of trees in tree1 and tree2 must match");
  }

  rspr_set_defaults();

  const bool COMPUTE_APPROX = approx[0];
  const bool COMPUTE_EXACT  = exact[0];

  // Label maps survive across tree pairs so that all pairs share a consistent
  // integer–label assignment (same as uspr.cpp).
  map<string, int> label_map;
  map<int, string> reverse_label_map;

  IntegerVector approx_dist(tree1.size(), NA_INTEGER);
  IntegerVector exact_dist (tree1.size(), NA_INTEGER);
  StringVector  maf_1      (tree1.size());
  StringVector  maf_2      (tree1.size());

  CoutSuppressor suppress;

  for (int i = 0; i < tree1.size(); i++) {
    const string tr1 = as<string>(tree1(i));
    const string tr2 = as<string>(tree2(i));

    Node *T1 = build_tree(tr1);
    Node *T2 = build_tree(tr2);

    T1->labels_to_numbers(&label_map, &reverse_label_map);
    T2->labels_to_numbers(&label_map, &reverse_label_map);

    // Preorder numbering is required by the clustering algorithm.
    // T2's preorder is also checked inside rSPR_branch_and_bound_simple_
    // clustering, but we call it here for the approx path too.
    T1->preorder_number();
    T1->edge_preorder_interval();
    T2->preorder_number();
    T2->edge_preorder_interval();

    if (COMPUTE_APPROX) {
      // rSPR_worse_3_approx is destructive and its RETURN VALUE is not the
      // distance.  The approx distance is F2.num_components() - 1.
      // sync_twins is called internally when sync = true (default).
      Forest F1 = Forest(T1);
      Forest F2 = Forest(T2);
      rSPR_worse_3_approx(&F1, &F2);
      approx_dist[i] = F2.num_components() - 1;
    }

    if (COMPUTE_EXACT) {
      Forest *MAF1 = NULL;
      Forest *MAF2 = NULL;
      exact_dist[i] = rSPR_branch_and_bound_simple_clustering(
        T1, T2, /*verbose=*/false,
        &label_map, &reverse_label_map,
        /*min_k=*/0, /*max_k=*/MAX_SPR,
        &MAF1, &MAF2);
      // Returns -1 if distance > MAX_SPR; NA_INTEGER is already the default.

      if (MAF1 != NULL) {
        MAF1->numbers_to_labels(&reverse_label_map);
        maf_1[i] = MAF1->str();
        delete MAF1;
      }
      if (MAF2 != NULL) {
        MAF2->numbers_to_labels(&reverse_label_map);
        maf_2[i] = MAF2->str();
        delete MAF2;
      }
    }

    T1->delete_tree();
    T2->delete_tree();

    Rcpp::checkUserInterrupt();
  }

  return List::create(
    Named("exact")  = exact_dist,
    Named("approx") = approx_dist,
    Named("maf_1")  = maf_1,
    Named("maf_2")  = maf_2
  );
}
