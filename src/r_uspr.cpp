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
int uspr_uspr (StringVector tree1,
               StringVector tree2,
               LogicalVector printMafs,
               LogicalVector countMafs,
               LogicalVector keepLabels,
               LogicalVector noOpt,
               LogicalVector protectB,
               LogicalVector tbrApprox,
               LogicalVector tbr,
               LogicalVector replug,
               LogicalVector uspr,
               LogicalVector approxEstimate,
               LogicalVector tbrEstimate,
               LogicalVector replugEstimate
               ) {
  /* protectB and *Estimate default to TRUE, all others to FALSE */
  PRINT_mAFS = printMafs[0];
  COUNT_mAFS = countMafs[0];
  KEEP_LABELS = keepLabels[0];
  DEFAULT_OPTIMIZATIONS = !noOpt[0];
  OPTIMIZE_PROTECT_B = protectB[0];
  COMPUTE_TBR_APPROX = tbrApprox[0];
  COMPUTE_TBR = tbr[0];
  COMPUTE_REPLUG = replug[0];
  COMPUTE_USPR = uspr[0];
  USE_TBR_APPROX_ESTIMATE = approxEstimate[0];
  USE_TBR_ESTIMATE = tbrEstimate[0];
  USE_REPLUG_ESTIMATE = replugEstimate[0];

  ALL_DISTANCES = !(COMPUTE_TBR_APPROX || COMPUTE_TBR || COMPUTE_REPLUG
                      || COMPUTE_USPR);


  if (DEFAULT_OPTIMIZATIONS == false) {
    OPTIMIZE_2B = false;
    OPTIMIZE_PROTECT_A = false;
    OPTIMIZE_PROTECT_B = false;
    OPTIMIZE_BRANCH_AND_BOUND = false;
    Rcout << "NO OPTIMIZATIONS" << endl;
  }

  if (ALL_DISTANCES) {
    COMPUTE_TBR_APPROX = true;
    COMPUTE_TBR = true;
    COMPUTE_REPLUG = true;
    COMPUTE_USPR = true;
  }

  // label maps to allow string labels
  map<string, int> label_map = map<string, int>();
  map<int, string> reverse_label_map = map<int, string>();

  // set random seed
  srand(unsigned(time(0)));

  // read input trees
  if (tree1.size() != tree2.size()) {
    throw length_error("Number of trees in tree1 and tree2 must match");
  }
  for (int i = 0; i < tree1.size(); i++) {
    // load into data structures
    uforest F1 = uforest(tree1[i], &label_map, &reverse_label_map);
    F1.normalize_order();
    uforest F2 = uforest(tree2[i], &label_map, &reverse_label_map);
    F2.normalize_order();
    Rcout << "T1: " << F1.str(false, &reverse_label_map) << endl;
    Rcout << "T2: " << F2.str(false, &reverse_label_map) << endl;
    // compute TBR distance
    if (COMPUTE_TBR_APPROX) {
      Rcout << "a_TBR: " << tbr_high_lower_bound(F1, F2) << " <= d_TBR <= " << tbr_low_upper_bound(F1, F2) << endl;
    }
    if (COMPUTE_TBR) {
      uforest *MAF1 = NULL;
      uforest *MAF2 = NULL;
      int distance = tbr_distance(F1, F2, false, &MAF1, &MAF2);
      Rcout << "d_TBR = " << distance << endl;
      if (MAF1 != NULL) {
        Rcout << "F1: " << MAF1->str(false, &reverse_label_map) << endl;
        delete MAF1;
      }
      if (MAF2 != NULL) {
        Rcout << "F2: " << MAF2->str(false, &reverse_label_map) << endl;
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

    if (COMPUTE_REPLUG) {
      uforest *MAF1 = NULL;
      uforest *MAF2 = NULL;
      int d_replug = replug_distance(F1, F2, false, &MAF1, &MAF2);
      Rcout << "d_R = " << d_replug << endl;
      if (MAF1 != NULL) {
        Rcout << "F1: " << MAF1->str(false, &reverse_label_map) << endl;
        delete MAF1;
      }
      if (MAF2 != NULL) {
        Rcout << "F2: " << MAF2->str(false, &reverse_label_map) << endl;
        delete MAF2;
      }
    }

    if (COMPUTE_USPR) {
      int d_uspr = uspr_distance(F1, F2);
      Rcout << "d_USPR = " << d_uspr << endl;
    }

  }
}
