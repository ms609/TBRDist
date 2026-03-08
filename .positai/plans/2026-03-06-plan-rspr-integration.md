# Plan: Integrate rspr as Rcpp submodule in TBRDist

## Background

**On alternatives:** rspr (Chris Whidden) is the state of the art for exact rooted SPR distance. Its O(2^k n) FPT algorithm with cluster decomposition is the best known. No better alternatives exist.

**What rspr provides over uspr:** rooted trees, 200+ leaves, distances >50. Complements TBRDist's existing unrooted-tree functionality.

## Key technical findings

- rspr is header-heavy (rspr.h is 7000+ lines) with all function bodies defined in headers, plus ~55 global variables (booleans/ints controlling algorithm options) also defined in headers. This is the same pattern as uspr.
- **Single translation unit (TU) approach** is required: a single `src/rspr.cpp` includes all rspr headers. This avoids multiple-definition errors from the header-defined globals, exactly as `src/uspr.cpp` does for uspr.
- rspr headers use `cout` extensively (452× in rspr.h alone). We suppress this with a `std::streambuf` swap inside the Rcpp wrapper — no modification to rspr source needed.
- `rSPR_branch_and_bound_simple_clustering(T1, T2, verbose=false, &lmap, &rlmap, 0, MAX_SPR, &F1_out, &F2_out)` is the main workhorse: exact distance with cluster reduction (the default and best algorithm). It requires `T1->preorder_number()` and `T1->edge_preorder_interval()` to be called first.
- `rSPR_worse_3_approx(&F1, &F2)` computes the 3-approximation. It is **destructive** — it modifies the Forest in place. The approx distance is `F2.num_components() - 1` (not the return value, which is 3× an internal counter).
- `Forest::str()` returns a space-separated string of Newick subtrees, e.g. `"(1,2) 3 (4,5) "` — this is the MAF representation.
- `Forest::numbers_to_labels(&reverse_label_map)` converts integer labels back to original strings.
- Trees are read from Newick strings via `build_tree(newick_string)` → `Node*`, then label-mapped with `T->labels_to_numbers(&label_map, &reverse_label_map)`.
- Memory: rspr uses raw pointers; Nodes are freed with `T->delete_tree()`, MAF Forests with `delete F`. Must ensure cleanup even on error/interrupt.
- rspr and uspr headers have **no naming conflicts** (rspr uses `Node`/`Forest`, uspr uses `unode`/`utree`/`uforest`). The two `.cpp` wrapper files are separate TUs.
- `src/rspr/rspr.cpp` (which has `main()`) lives in a subdirectory and is NOT compiled by R's build system.

## Files to create / modify

| File | Action |
|------|--------|
| `src/rspr.cpp` | **Create** — Rcpp wrapper |
| `R/rspr.R` | **Create** — R user-facing functions |
| `tests/testthat/test-rspr.R` | **Create** — regression tests |
| `DESCRIPTION` | **Modify** — Description, Authors, Copyright |
| `src/Makevars` | **Create if needed** — only if C++17 is required |

NAMESPACE and RcppExports.{R,cpp} are auto-generated.

---

## 1. `src/rspr.cpp` — Rcpp wrapper

```cpp
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#include <sstream>
#include <map>
#include <string>
#include <climits>

// rspr headers — all included in this one TU
#include "rspr/Node.h"
#include "rspr/Forest.h"
#include "rspr/ClusterForest.h"
#include "rspr/LCA.h"
#include "rspr/ClusterInstance.h"
#include "rspr/SiblingPair.h"
#include "rspr/UndoMachine.h"
#include "rspr/rspr.h"

// Reset all rspr globals to the same defaults as rspr.cpp's DEFAULT_OPTIMIZATIONS=true block.
// Must be called at the start of every Rcpp-exported function because these are
// process-level globals that persist across R-level calls within a session.
static void rspr_set_defaults() {
    BB = true;
    CLUSTER_TEST = true;
    PREFER_RHO = true;
    // Optimizations (from rspr.cpp main(), DEFAULT_OPTIMIZATIONS=true block)
    CUT_ALL_B                  = true;
    CUT_ONE_B                  = true;
    REVERSE_CUT_ONE_B          = true;
    REVERSE_CUT_ONE_B_3        = true;
    CUT_TWO_B                  = true;
    CUT_AC_SEPARATE_COMPONENTS = true;
    EDGE_PROTECTION            = true;
    EDGE_PROTECTION_TWO_B      = true;
    NEAR_PREORDER_SIBLING_PAIRS= true;
    PREORDER_SIBLING_PAIRS     = true;
    LEAF_REDUCTION             = true;
    LEAF_REDUCTION2            = true;
    PREFER_NONBRANCHING        = true;
    DEEPEST_PROTECTED_ORDER    = true;
    DEEPEST_ORDER              = true;
    APPROX_CUT_ONE_B           = true;
    APPROX_CUT_TWO_B           = true;
    APPROX_REVERSE_CUT_ONE_B   = true;
    CLUSTER_TUNE               = 30;
    // Non-default features off
    MULTIFURCATING             = false;
    ALL_MAFS                   = false;
    VERBOSE                    = false;
    MAX_SPR                    = 1000;
    CLUSTER_MAX_SPR            = MAX_SPR;
}

// [[Rcpp::export]]
List rspr_dist(const StringVector tree1,
               const StringVector tree2,
               const LogicalVector approx,
               const LogicalVector exact) {

    if (tree1.size() != tree2.size()) {
        throw std::length_error("Number of trees in tree1 and tree2 must match");
    }

    rspr_set_defaults();

    bool COMPUTE_APPROX = approx[0];
    bool COMPUTE_EXACT  = exact[0];

    map<string, int> label_map;
    map<int, string> reverse_label_map;

    IntegerVector approx_dist(tree1.size(), NA_INTEGER);
    IntegerVector exact_dist(tree1.size(),  NA_INTEGER);
    StringVector  maf_1(tree1.size());
    StringVector  maf_2(tree1.size());

    // Redirect cout to suppress rspr's internal progress output
    ostringstream discard;
    streambuf *old_cout = cout.rdbuf(discard.rdbuf());

    for (int i = 0; i < tree1.size(); i++) {
        string tr1 = as<string>(tree1(i));
        string tr2 = as<string>(tree2(i));

        Node *T1 = build_tree(tr1);
        Node *T2 = build_tree(tr2);
        T1->labels_to_numbers(&label_map, &reverse_label_map);
        T2->labels_to_numbers(&label_map, &reverse_label_map);

        // Preorder numbering required by the clustering algorithm
        T1->preorder_number();
        T1->edge_preorder_interval();
        T2->preorder_number();
        T2->edge_preorder_interval();

        if (COMPUTE_APPROX) {
            Forest F1 = Forest(T1);
            Forest F2 = Forest(T2);
            sync_twins(&F1, &F2);
            rSPR_worse_3_approx(&F1, &F2);          // destructive
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

    cout.rdbuf(old_cout);  // restore cout

    List ret = List::create(
        Named("exact")  = exact_dist,
        Named("approx") = approx_dist,
        Named("maf_1")  = maf_1,
        Named("maf_2")  = maf_2
    );
    return ret;
}
```

**Notes:**
- `cout.rdbuf()` is restored even if loop completes normally. For safety against exceptions/interrupts, a RAII guard struct would be cleaner — worth doing.
- The `label_map` accumulates across pairs (same as in uspr.cpp) which is fine: integers assigned to labels are consistent across all pairs in the batch.
- `rSPR_branch_and_bound_simple_clustering` returns -1 if the distance exceeds `MAX_SPR`; we return `NA_INTEGER` in that case (initialized as such).

---

## 2. `R/rspr.R` — user-facing wrapper

Single exported function `RSPRDist()`, modelled on `TBRDist()` / `USPRDist()`:

```r
#' Calculate rooted SPR distance between pairs of rooted trees
#'
#' @param tree1,tree2 Rooted trees of class `phylo`, or lists thereof.
#' @param allPairs Logical; compare all pairs or corresponding pairs.
#' @param checks Logical; validate input.
#' @param approx Logical; return 3-approximation only (faster). 
#' @param maf Logical; return maximum agreement forests.
#'
#' @return Integer vector of rSPR distances, or a list including MAFs when
#'   `maf = TRUE`.
#'
#' @references
#'   Whidden, Beiko & Zeh (2013) Fixed-Parameter Algorithms for Maximum
#'   Agreement Forests. _SIAM J. Comput._ 42:1431–1466.
#'   doi:10.1137/110845045
#'
#' @export
RSPRDist <- function(tree1, tree2 = NULL, allPairs = is.null(tree2),
                     checks = TRUE, approx = FALSE, maf = FALSE) {
    trees <- .PrepareTrees(tree1, tree2, allPairs = allPairs,
                           checks = checks, unroot = FALSE)
    ret <- rspr_dist(trees[[1]], trees[[2]],
                     approx = approx, exact = !approx)
    if (maf) {
        .DistReturn(ret[["exact"]], trees[[1]], trees[[2]], allPairs)
        # return full list including MAFs
        ...
    } else {
        dist <- if (approx) ret[["approx"]] else ret[["exact"]]
        .DistReturn(dist, trees[[1]], trees[[2]], allPairs)
    }
}
```

Key design points:
- `unroot = FALSE` passed to `.PrepareTrees()` since rspr requires rooted trees.
- If `checks = TRUE`, warn (not error) when an unrooted tree is passed, since the algorithm can still run.
- `approx = FALSE` by default (exact is the normal use case and rspr handles larger trees efficiently).
- When `maf = TRUE`, return a list analogous to `TBRDist(..., maf = TRUE)`.

---

## 3. `tests/testthat/test-rspr.R`

Use known results from rspr's own test suite:
- `trees2.txt`: `T1=((((1,2),(3,4)),((5,6),(7,8))),(((9,10),(11,12)),((13,14),(15,16))))`, `T2=(((7,8),...)` → exact drSPR = 4
- `cluster_test`: `T1=(((x,((b1,b3),b2)),y),(f,(a,c)))`, `T2=(((x,y),f),((a,((b1,b2),b3)),c))` → exact drSPR = 3
- Two identical trees → distance = 0
- Known approx results

---

## 4. `DESCRIPTION` changes

- Add rspr to Description and URL fields.
- Add a second `person()` entry for Chris Whidden with `role = "cph"` for rspr (he is already credited for uspr; may need a separate entry or note).
- Update Copyright line.

---

## Sequencing

1. Create `src/rspr.cpp`
2. Test compilation with `Rcpp::compileAttributes()` + `devtools::load_all()`
3. Create `R/rspr.R` + `tests/testthat/test-rspr.R`
4. Iterate on any compilation issues (likely: missing include, global collision, cout not fully suppressed)
5. Update `DESCRIPTION`
6. Run `R CMD check`

The biggest uncertainty is compilation: rspr is complex and the single-TU approach should work, but include order may matter (Forest.h and Node.h have a circular dependency that rspr resolves via include guards; we must follow the same include order as rspr.cpp's main does).
