#' Calculate rSPR distance between rooted trees
#'
#' Calculate rooted Subtree Prune-and-Regraft (rSPR) distances between pairs
#' of rooted binary trees using the algorithms of Whidden, Beiko & Zeh (2013).
#'
#' This function wraps the \pkg{rspr} C++ library by Chris Whidden. It handles
#' rooted trees and is substantially more efficient than
#' \code{\link{USPRDist}()} for large trees: it can handle pairs with 200+
#' leaves and distances > 50.
#'
#' Note that these distances are NP-hard to compute exactly; running time
#' scales as O(2^k n) where k is the rSPR distance and n is the number of
#' leaves. The built-in cluster decomposition (enabled by default) provides a
#' large practical speedup. Use \code{approx = TRUE} for a guaranteed
#' linear-time 3-approximation.
#'
#' Input trees must be **rooted**. An error is raised if any tree is unrooted.
#'
#' @param tree1,tree2 Trees of class \code{phylo}, or lists thereof.
#'   All trees must be rooted and bear identical sets of tip labels.
#' @param allPairs Logical; if \code{TRUE}, compare each tree in \code{tree1}
#'   with each tree in \code{tree2}; if \code{FALSE}, compare corresponding
#'   pairs.  Defaults to \code{TRUE} when \code{tree2} is not supplied.
#' @param checks Logical; validate tree labels and dimensions before
#'   computation.  Set \code{FALSE} at your peril—improper input is likely
#'   to crash R.
#' @param approx Logical; if \code{TRUE} return the linear-time
#'   3-approximation instead of the exact distance.
#' @param maf Logical; if \code{TRUE} return the maximum agreement forest
#'   alongside the exact distance (implies \code{approx = FALSE}).
#'
#' @return
#' By default, an integer vector of rSPR distances (one per tree pair), or a
#' \code{dist} object when \code{allPairs = TRUE} and \code{tree2 = NULL}.
#'
#' If \code{maf = TRUE}, a named list with elements:
#' \describe{
#'   \item{\code{exact}}{Integer vector of rSPR distances.}
#'   \item{\code{maf_1}, \code{maf_2}}{Character vectors giving the maximum
#'     agreement forest for each pair of trees, expressed as space-separated
#'     Newick components.}
#' }
#'
#' @examples
#' library(ape)
#' set.seed(1)
#' tree1 <- rtree(8)
#' tree2 <- rtree(8)
#' tree1$tip.label <- tree2$tip.label <- paste0("t", 1:8)
#'
#' RSPRDist(tree1, tree2)
#'
#' # All pairwise distances among a list of trees
#' trees <- c(tree1, tree2)
#' RSPRDist(trees)
#'
#' # Fast 3-approximation
#' RSPRDist(tree1, tree2, approx = TRUE)
#'
#' # With maximum agreement forest
#' RSPRDist(tree1, tree2, maf = TRUE)
#'
#' @references
#' Whidden C, Beiko RG, Zeh N (2013).
#' Fixed-Parameter Algorithms for Maximum Agreement Forests.
#' _SIAM Journal on Computing_ **42**:1431–1466.
#' \doi{10.1137/110845045}
#'
#' Whidden C, Zeh N, Beiko RG (2014).
#' Supertrees based on the subtree prune-and-regraft distance.
#' _Systematic Biology_ **63**:566–581.
#' \doi{10.1093/sysbio/syu023}
#'
#' @seealso \code{\link{USPRDist}()} for unrooted trees.
#' @importFrom ape is.rooted write.tree
#' @export
RSPRDist <- function(tree1, tree2 = NULL, allPairs = is.null(tree2),
                     checks = TRUE, approx = FALSE, maf = FALSE) {
  if (maf && approx) {
    warning("maf = TRUE requires exact computation; ignoring approx = TRUE")
    approx <- FALSE
  }

  if (checks) {
    .CheckRooted(tree1)
    if (!is.null(tree2)) .CheckRooted(tree2)
  }

  # write.tree produces rooted Newick with original tip labels, which rspr
  # parses correctly.  as.Newick (used by .PrepareTrees with keepLabels = FALSE)
  # may produce unrooted topology — so we always request keepLabels = TRUE here.
  treeLists <- .PrepareTrees(tree1, tree2, allPairs, checks,
                             unroot = FALSE, keepLabels = TRUE)

  ret <- rspr_dist(treeLists[[1]], treeLists[[2]],
                   approx = approx, exact = !approx)

  if (maf) {
    out <- list(exact = ret[["exact"]],
                maf_1 = ret[["maf_1"]],
                maf_2 = ret[["maf_2"]])
    .DistReturn(out, tree1, tree2, allPairs)
  } else {
    dist <- if (approx) ret[["approx"]] else ret[["exact"]]
    .DistReturn(dist, tree1, tree2, allPairs)
  }
}

#' Check that trees are rooted
#'
#' @param trees A \code{phylo} object or a list of \code{phylo} objects.
#' @return Called for its side-effect (stops with an error if any tree is
#'   unrooted).
#' @keywords internal
.CheckRooted <- function(trees) {
  if (inherits(trees, "phylo")) trees <- list(trees)
  rooted <- vapply(trees, is.rooted, logical(1))
  if (!all(rooted)) {
    stop("RSPRDist() requires rooted trees. ",
         "Use ape::root() or phytools::midpoint.root() to root your trees, ",
         "or USPRDist() if you want unrooted SPR distances.")
  }
}
