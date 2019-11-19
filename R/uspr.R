#' @export
#' @importFrom ape write.tree
USPRDist <- function (tree1, tree2) {
  r_uspr(write.tree(tree1), write.tree(tree2),
         printMafs = FALSE, countMafs = FALSE, keepLabels = TRUE,
         opt = TRUE, protectB = TRUE,
         tbrApprox = FALSE, tbr = FALSE, replug = FALSE,
         uspr = TRUE,
         approxEstimate = TRUE, tbrEstimate = TRUE, replugEstimate = TRUE)
}
