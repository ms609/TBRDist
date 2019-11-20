#' Calculate Tree Rearrangement Distances
#'
#' Calculate SPR, TBR and Replug distances on unrooted trees.
#'
#'
#' Note that these distances are NP-hard to compute and most of the algorithms
#' used in this software scale exponentially with the distance computed.
#' This version of uspr is aimed at trees with up to 50 leaves and uSPR
#' distances up to 14.
#'
#' If you are interested in comparing rooted trees in terms of SPR operations,
#' you should use [rspr](https://github.com/cwhidden/rspr) instead. rspr is also
#' much more efficient and can easily handle pairs of binary rooted trees with
#' 200+ leaves and distances > 50.
#' rspr is not yet incorporated in this R package; please
#' [contact the maintainer](https://github.com/ms609/uspr/issues/2)
#' if this would be useful to you.
#'
#'
#' @param tree1,tree2 Trees of class `phylo`, or lists thereof.
#' @param checks Logical specifying whether to check that trees are properly
#' formatted and labelled.  Specify `FALSE` at your peril, as improper
#' input is likely to crash R.
#' @param allPairs Logical; if `TRUE`, compare each tree in `tree1` with each
#' tree in `tree2`; if `FALSE`, compare each tree in `tree1` only with the
#' tree at the corresponding index in `tree2`.  If `tree2` is not specified,
#' each tree in `tree1` will be compared with each other tree in `tree1`.
#'
#' @param useTbrApproxEstimate,useTbrEstimate,useReplugEstimate Logical specifying
#' whether to use approximate TBR distance, TBR distance or Replug distance to
#' help estimate the SPR distance.
#'
#' @return `USPRDist` returns a vector of SPR distances between each pair of
#' unrooted trees.
#'
#' @examples
#'   tree1 <- TreeTools::BalancedTree(9)
#'   tree2 <- TreeTools::PectinateTree(9)
#'
#'   # SPR distance
#'   USPRDist(tree1, tree2)
#'
#'
#' @name TreeRearrangementDistances
#' @author
#' Algorithms implemented by Chris Whidden (<cwhidden@fredhutch.org>)
#'
#' R wrappers by Martin R. Smith (<martin.smith@durham.ac.uk>)
#'
#' @references
#' If you use these functions in your research, please cite:
#'
#' * Chris Whidden and Frederick A. Matsen IV. Calculating the Unrooted
#' Subtree-Prune-and-Regraft Distance.
#' arXiv:[1511.07529](http://arxiv.org/abs/1511.07529).
#'
#' @export
USPRDist <- function (tree1, tree2 = NULL, allPairs = is.null(tree2),
                      checks = TRUE,
                      useTbrApproxEstimate = TRUE,
                      useTbrEstimate = TRUE,
                      useReplugEstimate = TRUE) {
  treeLists <- .PrepareTrees(tree1, tree2, allPairs, checks)
  ret <- uspr_dist(treeLists[[1]], treeLists[[2]], keepLabels = FALSE,
                   useTbrApproxEstimate = useTbrApproxEstimate,
                   useTbrEstimate = useTbrEstimate,
                   useReplugEstimate = useReplugEstimate)
  .DistReturn(ret, tree1, tree2, allPairs)
}

#' @rdname TreeRearrangementDistances
#' @param maf Logical specifying whether to report a maximum agreement forest
#' corresponding to the optimal score.
#'
#' @return `ReplugDist` returns a vector of Replug distances between each pair
#' of trees, or (if `maf = TRUE`) a named list whose second and third elements
#' list a vector of maximum agreement forests for each pair of trees.
#'
#' @examples
#'   # Replug distance
#'   ReplugDist(tree1, tree2)
#'   ReplugDist(tree1, tree2, maf = TRUE)
#'
#' @export
ReplugDist <- function (tree1, tree2 = NULL, allPairs = is.null(tree2),
                        checks = TRUE, maf = FALSE) {
  treeLists <- .PrepareTrees(tree1, tree2, allPairs, checks)
  ret <- replug_dist(treeLists[[1]], treeLists[[2]], keepLabels = FALSE)
  if (maf) {
    names(ret) <- c('replug_dist', 'maf_1', 'maf_2')
    .DistReturn(ret, tree1, tree2, allPairs)
  } else {
    .DistReturn(ret[[1]], tree1, tree2, allPairs)
  }
}

#' @rdname TreeRearrangementDistances
#' @param exact Logical specifying whether to calculate the exact TBR distance.
#' @param approximate Logical specifying whether to calculate the approximate
#' TBR distance.  By default, is set to the opposite of `exact`; either
#' `approximate` or `exact` should usually be set to `TRUE` if a distance is
#' required.
#' @param countMafs Logical specifying whether to count the number of Maximum
#' Agreement Forests found.
#' @param printMafs Logical specifying whether to print Maximum Agreement
#' Forests to stdout whilst counting.
#' Use [`capture.output`]`(TBRDist(tree1, tree2, printMafs = TRUE))` to access
#' these in R.
#' @param optimize Logical specifying whether to use the default optimizations.
#' @param protectB Logical specifying whether to use the PROTECT_B optimization.
#' Overrides value inherited from `optimize`.
#' @return `TBRDist` returns a named list, each element of which bears a vector
#' corresponding to the requested value for each tree pair.  If only the exact
#' value is requested (`exact = TRUE`), an unnamed vector of distances is
#' returned.
#'
#' @examples
#'   # TBR distance between two trees
#'   TBRDist(tree1, tree2, exact = TRUE)
#'
#'   # Compare a list against one tree, using approximate distances
#'   TBRDist(list(tree1, tree2), tree2, exact = FALSE)
#'
#'   # Compare all pairs in two lists
#'   TBRDist(list(tree1, tree2), list(tree1, tree2, tree2), allPairs = TRUE,
#'           exact = FALSE)
#'
#'   # Compare each tree in a list against each other
#'   TBRDist(list(one = tree1, two = tree2, twoAgain = tree2))
#'
#'   # Compare each pair in two lists
#'   TBRDist(list(tree1, tree2, tree2),
#'           list(tree2, tree1, tree2),
#'           exact = TRUE, approximate = TRUE, countMafs = TRUE)
#'
#'   # Capture maximum agreement forests
#'   mafs <- capture.output(TBRDist(tree1, tree2, approximate = FALSE,
#'                           printMafs = TRUE))
#'   head(mafs)
#'
#' @export
TBRDist <- function (tree1, tree2 = NULL, allPairs = is.null(tree2),
                     checks = TRUE,
                     exact = FALSE, approximate = !exact,
                     maf = FALSE, countMafs = FALSE, printMafs = FALSE,
                     optimize = TRUE, protectB = TRUE) {
  if (!exact && !approximate && !countMafs && !printMafs) {
    message("Nothing to do in TBRDist.")
  }
  treeLists <- .PrepareTrees(tree1, tree2, allPairs, checks)

  whichRets <- c(exact, rep(approximate, 2L), countMafs,
                 rep(exact && maf, 2L))

  ret <- tbr_dist(treeLists[[1]], treeLists[[2]],
                  printMafs = printMafs, countMafs = countMafs,
                  keepLabels = FALSE,
                  optimize = optimize, protectB = protectB,
                  exact = exact, approximate = approximate)[whichRets]
  if (!any(whichRets)) {
    invisible()
  } else if (exact && sum(whichRets) == 1) {
    .DistReturn(ret[[1]], tree1, tree2, allPairs)
  } else {
    names(ret) <- c('tbr_exact', 'tbr_min', 'tbr_max', 'n_maf', 'maf_1', 'maf_2')[whichRets]
    .DistReturn(ret, tree1, tree2, allPairs)
  }
}

#' @keywords internal
#' @importFrom ape write.tree
#' @export
.PrepareTrees <- function (tree1, tree2, allPairs = FALSE, checks = TRUE) {
  if (checks) {

    if (class(tree1) == 'phylo') tree1 <- list(tree1)
    if (class(tree2) == 'phylo') tree2 <- list(tree2)

    if (allPairs) {
      if (is.null(tree2)) {
        nTree <- length(tree1)
        selector <- matrix(seq_len(nTree), nTree, nTree)
        tree2 <- tree1[t(selector)[lower.tri(selector)]]
        tree1 <- tree1[selector[lower.tri(selector)]]
      } else {
        nTree1 <- length(tree1)
        tree1 <- rep(tree1, each = length(tree2))
        tree2 <- rep(tree2, nTree1)
      }
    } else {
      if (length(tree1) != length(tree2)) {
        if (length(tree1) == 1L) {
          tree1 <- rep(tree1, length(tree2))
        } else if (length(tree2) == 1L) {
          tree2 <- rep(tree2, length(tree1))
        } else {
          stop("if allPairs = FALSE, tree1 and tree2 must contain the same number of trees, or a single tree.")
        }
      }
    }

    vapply(seq_along(tree1), function (i) .CatchBadPair(i, tree1[[i]], tree2[[i]]),
           integer(0))
  }
  class(tree1) <- class(tree2) <- 'multiPhylo'
  list(write.tree(tree1), write.tree(tree2))
}

#' @keywords internal
#' @export
.CatchBadPair <- function (i, tree1, tree2) {
  lab1 <- tree1$tip.label
  lab2 <- tree2$tip.label
  if (length(lab1) != length(lab2)) {
    stop("Problem with tree pair " , i, ": Trees must be the same size")
  }
  if (length(setdiff(tree1$tip.label, tree2$tip.label)) > 0){
    stop("Problem with tree pair " , i, ": Trees must bear identical labels")
  }
  integer(0)
}

#' @keywords internal
#' @export
.DistReturn <- function (ret, tree1, tree2, allPairs) {
  .ReturnMatrix <- function (dat) {
    matrix(dat, length(tree1), length(tree2), byrow = TRUE,
           dimnames = list(names(tree1), names(tree2)))
  }
  .ReturnDist <- function (dat, nTree) {
    if (mode(dat) == 'numeric') {
      ret <- structure(dat, Size = nTree, Diag = FALSE, Upper = FALSE,
                       Labels = names(tree1), class = 'dist')
    } else {
      ret <- matrix(NA, nTree, nTree)
      ret[lower.tri(ret)] <- dat
      ret[upper.tri(ret)] <- t(dat)
    }

    # Return:
    ret
  }
  if (allPairs) {
    if (is.null(tree2)) {
      if (mode(ret) == 'list') {
        lapply(ret, .ReturnDist, length(tree1))
      } else {
        .ReturnDist(ret, length(tree1))
      }
    } else {
      names1 <- names(tree1)
      names2 <- names(tree2)
      if (mode(ret) == 'list') {
        lapply(ret, .ReturnMatrix)
      } else {
        .ReturnMatrix(ret)
      }
    }
  } else {
    ret
  }
}
