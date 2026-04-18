# R/plot_critical_edges.R

#' Visualize Central Graphs and Critical Edges
#'
#' Produces a multi-panel figure summarising the output of
#' \code{\link{identify_critical_links}}: one panel per population showing the
#' weighted central graph, plus a final panel that highlights the critical
#' edges on a chosen reference central graph. This is the canonical
#' visualization of the \pkg{BrainNetTest} workflow and replaces the manual
#' \pkg{igraph} plumbing required to build the same figure.
#'
#' @param populations A named \code{list} of populations of adjacency matrices,
#'   as accepted by \code{\link{identify_critical_links}}. Used to compute the
#'   per-population central graphs unless \code{central_graphs} is supplied.
#' @param critical_links The object returned by
#'   \code{\link{identify_critical_links}}. Its \code{critical_edges} component
#'   provides the edges to highlight in the final panel.
#' @param communities Optional integer (or factor) vector of length
#'   \code{n_nodes} giving the community membership of every node. Used to
#'   color vertices consistently across panels. Default \code{NULL} (all
#'   vertices share the default \pkg{igraph} color).
#' @param layout Optional layout for the graphs. Either an
#'   \code{n_nodes x 2} numeric matrix of vertex coordinates (recommended, so
#'   that all panels share the same node positions) or an \pkg{igraph} layout
#'   function such as \code{\link[igraph]{layout_with_fr}}. Default
#'   \code{\link[igraph]{layout_in_circle}}.
#' @param central_graphs Optional list of pre-computed central graphs (one per
#'   population). If \code{NULL}, they are computed via
#'   \code{\link{compute_central_graph}}.
#' @param reference Either an integer index or a name selecting which
#'   population's central graph is used as the background of the
#'   critical-edges panel. Default \code{1L} (the first population).
#' @param threshold Numeric in \code{[0, 1]}. The reference central graph is
#'   binarised by \code{> threshold} for the critical-edges panel. Default
#'   \code{0.5}.
#' @param edge_scale Numeric multiplier applied to the central-graph edge
#'   weights to set \code{edge.width} in the per-population panels. Default
#'   \code{6}.
#' @param critical_color,critical_width Color and line width used to draw
#'   critical edges. Defaults \code{"red"} and \code{3}.
#' @param background_color,background_width Color and line width used to draw
#'   non-critical edges in the critical panel. Defaults \code{"grey80"} and
#'   \code{1}.
#' @param vertex_size Vertex size passed to \code{\link[igraph]{plot.igraph}}.
#'   Default \code{12}.
#' @param vertex_label Logical. If \code{TRUE} (default) vertex indices are
#'   shown; otherwise vertex labels are suppressed.
#' @param panel_titles Optional character vector of length
#'   \code{length(populations) + 1} with custom panel titles. Default uses
#'   \code{paste(names(populations), "Central Graph")} followed by
#'   \code{"Critical edges"}.
#' @param mfrow Optional integer vector of length 2 giving the
#'   \code{c(nrow, ncol)} arrangement passed to \code{\link[graphics]{par}}.
#'   By default a near-square grid is chosen automatically (e.g.\ a
#'   \code{2 x 2} grid for two populations, matching the standard
#'   Control / Patient / Critical layout).
#' @param ... Additional arguments forwarded to
#'   \code{\link[igraph]{plot.igraph}}. Note that \code{edge.width},
#'   \code{vertex.color}, \code{vertex.size}, \code{vertex.label} and
#'   \code{layout} are managed by this function and should not be passed via
#'   \code{...}.
#'
#' @return Invisibly returns \code{NULL}; called for its side effect of
#'   producing the multi-panel plot.
#'
#' @details The first \code{length(populations)} panels show the weighted
#'   central graph of each population, with edge widths proportional to the
#'   central-graph weights. The final panel binarises the chosen reference
#'   central graph at \code{threshold} and overlays the critical edges
#'   returned by \code{identify_critical_links()} in
#'   \code{critical_color}. Critical edges that are absent from the binarised
#'   reference graph are added so that they remain visible.
#'
#'   The graphics state (\code{par(mfrow, mar, oma)}) is restored on exit.
#'
#' @seealso \code{\link{identify_critical_links}},
#'   \code{\link{compute_central_graph}},
#'   \code{\link{get_critical_nodes}}.
#'
#' @export
#' @importFrom graphics par
#' @examples
#' \donttest{
#' set.seed(123)
#' community_sizes <- c(4, 2, 3, 3, 5)
#' control <- generate_category_graphs(
#'   n_graphs = 50, n_nodes = 17, n_communities = 5,
#'   community_sizes = community_sizes,
#'   base_intra_prob = rep(0.70, 5), base_inter_prob = 0.05,
#'   intra_prob_variation = 0.02, inter_prob_variation = 0.01)
#' patient <- generate_category_graphs(
#'   n_graphs = 50, n_nodes = 17, n_communities = 5,
#'   community_sizes = community_sizes,
#'   base_intra_prob = c(0.40, 0.70, 0.70, 0.70, 0.70),
#'   base_inter_prob = 0.05,
#'   intra_prob_variation = 0.02, inter_prob_variation = 0.01)
#' populations <- list(Control = control, Patient = patient)
#' result <- identify_critical_links(populations, alpha = 0.05,
#'   method = "fisher", n_permutations = 200, seed = 1)
#' communities <- rep(seq_along(community_sizes), times = community_sizes)
#' plot_critical_edges(populations, result, communities = communities)
#' }
plot_critical_edges <- function(populations,
                                critical_links,
                                communities = NULL,
                                layout = NULL,
                                central_graphs = NULL,
                                reference = 1L,
                                threshold = 0.5,
                                edge_scale = 6,
                                critical_color = "red",
                                critical_width = 3,
                                background_color = "grey80",
                                background_width = 1,
                                vertex_size = 12,
                                vertex_label = TRUE,
                                panel_titles = NULL,
                                mfrow = NULL,
                                ...) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required. ",
         "Install it with install.packages('igraph').")
  }
  if (!is.list(populations) || length(populations) < 1L) {
    stop("'populations' must be a non-empty list of populations.")
  }
  if (is.null(names(populations))) {
    names(populations) <- paste0("Group", seq_along(populations))
  }
  if (!is.list(critical_links) ||
      !("critical_edges" %in% names(critical_links))) {
    stop("'critical_links' must be the output of identify_critical_links().")
  }

  ## Central graphs ---------------------------------------------------------
  if (is.null(central_graphs)) {
    central_graphs <- lapply(populations, compute_central_graph)
  }
  if (length(central_graphs) != length(populations)) {
    stop("'central_graphs' must have the same length as 'populations'.")
  }
  n_nodes <- nrow(central_graphs[[1L]])
  if (any(vapply(central_graphs, nrow, integer(1L)) != n_nodes)) {
    stop("All central graphs must have the same number of nodes.")
  }

  ## Reference index --------------------------------------------------------
  if (is.character(reference)) {
    ref_idx <- match(reference, names(populations))
    if (is.na(ref_idx)) {
      stop("'reference' name not found in names(populations).")
    }
  } else {
    ref_idx <- as.integer(reference)
    if (length(ref_idx) != 1L || is.na(ref_idx) ||
        ref_idx < 1L || ref_idx > length(populations)) {
      stop("'reference' index out of range.")
    }
  }

  ## Communities ------------------------------------------------------------
  if (!is.null(communities) && length(communities) != n_nodes) {
    stop("'communities' must have length equal to the number of nodes.")
  }

  ## Layout -----------------------------------------------------------------
  if (is.null(layout)) {
    layout <- igraph::layout_in_circle
  }

  ## Titles -----------------------------------------------------------------
  K <- length(populations)
  if (is.null(panel_titles)) {
    panel_titles <- c(paste(names(populations), "Central Graph"),
                      "Critical edges")
  } else if (length(panel_titles) != K + 1L) {
    stop("'panel_titles' must have length equal to length(populations) + 1.")
  }

  ## Layout grid ------------------------------------------------------------
  if (is.null(mfrow)) {
    n_panels <- K + 1L
    nc <- ceiling(sqrt(n_panels))
    nr <- ceiling(n_panels / nc)
    mfrow <- c(nr, nc)
  }

  old_par <- graphics::par(mfrow = mfrow,
                           mar = c(1, 1, 2, 1),
                           oma = c(0, 0, 0, 0))
  on.exit(graphics::par(old_par), add = TRUE)

  set_vertex_attrs <- function(g) {
    if (!is.null(communities)) {
      igraph::V(g)$color <- communities
    }
    igraph::V(g)$size <- vertex_size
    if (!isTRUE(vertex_label)) {
      igraph::V(g)$label <- NA
    }
    g
  }

  ## Per-population weighted central-graph panels --------------------------
  for (k in seq_len(K)) {
    M <- central_graphs[[k]]
    g <- igraph::graph_from_adjacency_matrix(M, mode = "undirected",
                                             weighted = TRUE, diag = FALSE)
    g <- set_vertex_attrs(g)
    w <- igraph::E(g)$weight
    if (is.null(w)) w <- 1
    plot(g, layout = layout, main = panel_titles[k],
         edge.width = w * edge_scale, ...)
  }

  ## Critical-edges panel ---------------------------------------------------
  ref_central <- central_graphs[[ref_idx]]
  bin <- ref_central > threshold
  diag(bin) <- FALSE
  g <- igraph::graph_from_adjacency_matrix(bin, mode = "undirected",
                                           diag = FALSE)
  g <- set_vertex_attrs(g)
  igraph::E(g)$color <- background_color
  igraph::E(g)$width <- background_width

  ce <- critical_links$critical_edges
  if (!is.null(ce) && nrow(ce) > 0L) {
    el <- igraph::as_edgelist(g)
    for (i in seq_len(nrow(ce))) {
      a <- ce$node1[i]
      b <- ce$node2[i]
      idx <- which((el[, 1L] == a & el[, 2L] == b) |
                   (el[, 1L] == b & el[, 2L] == a))
      if (length(idx) > 0L) {
        igraph::E(g)[idx]$color <- critical_color
        igraph::E(g)[idx]$width <- critical_width
      } else {
        ## Critical edge missing from the thresholded reference -- add it
        ## so that it appears in the figure.
        g <- igraph::add_edges(
          g, c(a, b),
          attr = list(color = critical_color, width = critical_width))
      }
    }
  }

  plot(g, layout = layout, main = panel_titles[K + 1L], ...)

  invisible(NULL)
}
