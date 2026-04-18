# R/get_critical_nodes.R

#' Extract Critical Nodes from Critical Edge Results
#'
#' Given the output of \code{\link{identify_critical_links}}, this function
#' identifies the unique nodes involved in the critical edges and summarizes
#' their participation. Each node is reported with the number of critical edges
#' it participates in (its \emph{critical degree}). If node labels are supplied,
#' they are included in the output.
#'
#' @param result The list returned by \code{\link{identify_critical_links}}.
#'   Must contain a \code{critical_edges} component (a data frame with columns
#'   \code{node1} and \code{node2}, or \code{NULL}).
#' @param node_labels An optional character vector of length \eqn{p} (the number
#'   of nodes in the network), where \code{node_labels[i]} is the label for
#'   node \code{i}. If \code{NULL} (the default), no labels are included.
#'
#' @return A data frame with the following columns, ordered by decreasing
#'   \code{critical_degree}:
#'   \describe{
#'     \item{node}{Integer node index.}
#'     \item{label}{Character label for the node (only present when
#'       \code{node_labels} is not \code{NULL}).}
#'     \item{critical_degree}{Number of critical edges incident on this node.}
#'   }
#'   If no critical edges were found (\code{result$critical_edges} is
#'   \code{NULL} or has zero rows), an empty data frame with the same columns
#'   is returned.
#'
#' @export
#'
#' @examples
#' # Generate two synthetic populations with different community structure
#' control <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
#'   n_communities = 2, base_intra_prob = 0.8, base_inter_prob = 0.2, seed = 1)
#' patient <- generate_category_graphs(n_graphs = 15, n_nodes = 10,
#'   n_communities = 2, base_intra_prob = 0.5, base_inter_prob = 0.5, seed = 2)
#' populations <- list(Control = control, Patient = patient)
#'
#' # Identify critical edges
#' result <- identify_critical_links(populations, alpha = 0.05,
#'   method = "fisher", n_permutations = 200, seed = 42)
#'
#' # Extract critical nodes (no labels)
#' get_critical_nodes(result)
#'
#' # With labels
#' labels <- paste0("Region_", 1:10)
#' get_critical_nodes(result, node_labels = labels)
get_critical_nodes <- function(result, node_labels = NULL) {

  ## --- Input validation ---
  if (!is.list(result) || !"critical_edges" %in% names(result))
    stop("`result` must be the list returned by identify_critical_links().")

  if (!is.null(node_labels) && !is.character(node_labels))
    stop("`node_labels` must be a character vector or NULL.")

  ## --- Handle empty results ---
  ce <- result$critical_edges
  if (is.null(ce) || nrow(ce) == 0L) {
    out <- data.frame(node = integer(0), critical_degree = integer(0))
    if (!is.null(node_labels))
      out$label <- character(0)
    return(out)
  }

  ## --- Validate node_labels length against the data ---
  max_node <- max(ce$node1, ce$node2)
  if (!is.null(node_labels) && length(node_labels) < max_node)
    stop(sprintf(
      "`node_labels` has length %d but critical edges reference node %d.",
      length(node_labels), max_node))

  ## --- Compute critical degree per node ---
  all_nodes <- c(ce$node1, ce$node2)
  freq <- table(all_nodes)
  node_ids <- as.integer(names(freq))
  degrees  <- as.integer(freq)

  out <- data.frame(node = node_ids, critical_degree = degrees)

  if (!is.null(node_labels))
    out$label <- node_labels[out$node]

  ## --- Sort by critical_degree descending ---
  out <- out[order(out$critical_degree, decreasing = TRUE), ]
  rownames(out) <- NULL

  ## --- Reorder columns so label comes before degree if present ---
  if (!is.null(node_labels))
    out <- out[, c("node", "label", "critical_degree")]

  out
}
