#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <memory>
#include <ctime>
#include <cstdlib>
#include "uspr/utree.h"
#include "uspr/unode.h"
#include "uspr/uforest.h"
#include "uspr/tbr.h"
#include "uspr/uspr.h"

// [[Rcpp::export]]
IntegerVector uspr_dist (StringVector tree1,
                         StringVector tree2,
                         LogicalVector keepLabels) {
  /* opt, protectB and *Estimate default to TRUE, all others to FALSE */
  KEEP_LABELS = keepLabels[0];

  // label maps to allow string labels
  map<string, int> label_map = map<string, int>();
  map<int, string> reverse_label_map = map<int, string>();

  // read input trees
  if (tree1.size() != tree2.size()) {
    throw length_error("Number of trees in tree1 and tree2 must match");
  }
  IntegerVector ret(tree1.size());
  for (int i = 0; i < tree1.size(); i++) {
    // load into data structures
    string tr1 = as<string>(tree1(i));
    string tr2 = as<string>(tree2(i));
    uforest F1 = uforest(tr1, &label_map, &reverse_label_map);
    F1.normalize_order();
    uforest F2 = uforest(tr2, &label_map, &reverse_label_map);
    F2.normalize_order();
    ret(i) = uspr_distance(F1, F2);
  }
  return (ret);
}

// [[Rcpp::export]]
List tbr_dist (StringVector tree1,
                        StringVector tree2,
                        LogicalVector printMafs,
                        LogicalVector countMafs,
                        LogicalVector keepLabels,
                        LogicalVector opt,
                        LogicalVector protectB,
                        LogicalVector tbrApprox,
                        LogicalVector tbr,
                        LogicalVector approxEstimate,
                        LogicalVector tbrEstimate) {
  /* opt, protectB and *Estimate default to TRUE, all others to FALSE */
  bool PRINT_mAFS = printMafs[0];
  bool COUNT_mAFS = countMafs[0];
  KEEP_LABELS = keepLabels[0];
  bool DEFAULT_OPTIMIZATIONS = opt[0];
  OPTIMIZE_PROTECT_B = protectB[0];
  bool COMPUTE_TBR_APPROX = tbrApprox[0];
  bool COMPUTE_TBR = tbr[0];

  USE_TBR_APPROX_ESTIMATE = approxEstimate[0];
  USE_TBR_ESTIMATE = tbrEstimate[0];

  if (DEFAULT_OPTIMIZATIONS == false) {
    OPTIMIZE_2B = false;
    OPTIMIZE_PROTECT_A = false;
    OPTIMIZE_PROTECT_B = false;
    OPTIMIZE_BRANCH_AND_BOUND = false;
  }

  // label maps to allow string labels
  map<string, int> label_map = map<string, int>();
  map<int, string> reverse_label_map = map<int, string>();

  // read input trees
  if (tree1.size() != tree2.size()) {
    throw length_error("Number of trees in tree1 and tree2 must match");
  }
  IntegerVector tbr_exact(tree1.size()),
    tbr_above(tree1.size()),
    tbr_below(tree1.size());
  StringVector maf_1(tree1.size()),
    maf_2(tree1.size());
  for (int i = 0; i < tree1.size(); i++) {
    // load into data structures
    string tr1 = as<string>(tree1(i));
    string tr2 = as<string>(tree2(i));
    uforest F1 = uforest(tr1, &label_map, &reverse_label_map);
    F1.normalize_order();
    uforest F2 = uforest(tr2, &label_map, &reverse_label_map);
    F2.normalize_order();
    // compute TBR distance
    if (COMPUTE_TBR_APPROX) {
      tbr_above(i) = tbr_high_lower_bound(F1, F2);
      tbr_below(i) = tbr_low_upper_bound(F1, F2);
    }
    if (COMPUTE_TBR) {
      uforest *MAF1 = NULL;
      uforest *MAF2 = NULL;
      tbr_exact[i] = tbr_distance(F1, F2, false, &MAF1, &MAF2);

      if (MAF1 != NULL) {
        maf_1(i) = MAF1->str(false, &reverse_label_map);
        delete MAF1;
      }
      if (MAF2 != NULL) {
        maf_2(i) = MAF2->str(false, &reverse_label_map);
        delete MAF2;
      }
    }
    int count;
    if (PRINT_mAFS) {
      count = tbr_print_mAFs(F1, F2);
      Rcout << count << " mAFs" << endl;
    }
    else if (COUNT_mAFS) {
      count = tbr_count_mAFs(F1, F2);
      Rcout << count << " mAFs" << endl;
    }
  }
  List ret = List::create(tbr_exact, tbr_above, tbr_below, maf_1, maf_2);
  return (ret);
}

// [[Rcpp::export]]
IntegerVector replug_dist (StringVector tree1,
                           StringVector tree2,
                           LogicalVector keepLabels,
                           LogicalVector replugEstimate) {
  /* opt, replugEstimate defaults to TRUE */

  KEEP_LABELS = keepLabels[0];

  // label maps to allow string labels
  map<string, int> label_map = map<string, int>();
  map<int, string> reverse_label_map = map<int, string>();

  // read input trees
  if (tree1.size() != tree2.size()) {
    throw length_error("Number of trees in tree1 and tree2 must match");
  }
  IntegerVector replug(tree1.size());
  for (int i = 0; i < tree1.size(); i++) {
    // load into data structures
    string tr1 = as<string>(tree1(i));
    string tr2 = as<string>(tree2(i));
    uforest F1 = uforest(tr1, &label_map, &reverse_label_map);
    F1.normalize_order();
    uforest F2 = uforest(tr2, &label_map, &reverse_label_map);
    F2.normalize_order();

    uforest *MAF1 = NULL;
    uforest *MAF2 = NULL;
    replug(i) = replug_distance(F1, F2, false, &MAF1, &MAF2);

    if (MAF1 != NULL) {
      Rcout << "F1: " << MAF1->str(false, &reverse_label_map) << endl;
      delete MAF1;
    }
    if (MAF2 != NULL) {
      Rcout << "F2: " << MAF2->str(false, &reverse_label_map) << endl;
      delete MAF2;
    }

  }
  return (replug);
}
