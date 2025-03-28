#----- Prepare Data -----#
# Test 19: SingleEdgeDiff_3x3
# Base adjacency matrix
base_mat <- matrix(c(
  0, 1, 1,
  1, 0, 0,
  1, 0, 0
), nrow = 3, ncol = 3, byrow = TRUE)

# Population A: Force edge (2,3) = 1
A1 <- base_mat
A1[2, 3] <- 1
A1[3, 2] <- 1  # Keep symmetry
A2 <- A1  # Identical second graph

# Population B: Keep that edge at 0
B1 <- base_mat
B2 <- base_mat

# Create populations
populations <- list(A = list(A1, A2), B = list(B1, B2))

#--------------------#
#--------------------#
#-------Different Data-------------#
#--------------------#
#--------------------#
# Generate distinctly different graph populations
A <- generate_category_graphs(
  n_graphs = 5,
  n_nodes = 5,
  n_communities = 1,
  base_intra_prob = 1,
  base_inter_prob = 0.0,
  intra_prob_variation = 0.0,
  inter_prob_variation = 0.0,
  seed = 42
)

B <- generate_category_graphs(
  n_graphs = 5,
  n_nodes = 5,
  n_communities = 1,
  base_intra_prob = 0.1,
  base_inter_prob = 0.0,
  intra_prob_variation = 0.0,
  inter_prob_variation = 0.0,
  seed = 2
)

populations <- list(A = A, B = B)
#--------------------#
#--------------------#
#--------------------#
#--------------------#

#----- execute statistic test -----#


if (!is.list(populations) || length(populations) == 0) {
  stop("populations must be a non-empty list.")
}

m <- length(populations)  # Number of populations
n_i <- sapply(populations, length)  # Number of graphs in each population
n <- sum(n_i)  # Total number of graphs

if (any(n_i < 2)) {
  stop("Each population must have at least two graphs.")
}

# Compute central graphs for each population
# Since all graphs are the same central graph should be equal
central_graphs <- lapply(populations, compute_central_graph)

# Compute average distance of each population to its central graph
avg_d_Gi_Mi <- numeric(m)
names(avg_d_Gi_Mi) <- names(populations)

for (i in seq_along(populations)) {
  distances <- sapply(populations[[i]], compute_distance, M = central_graphs[[i]])
  avg_d_Gi_Mi[i] <- mean(distances)
}

# Compute average distance of all graphs to each central graph
avg_d_G_Mi <- numeric(m)
names(avg_d_G_Mi) <- names(populations)

for (i in seq_along(central_graphs)) {
  all_distances <- numeric(n)
  idx <- 1
  for (j in seq_along(populations)) {
    distances <- sapply(populations[[j]], compute_distance, M = central_graphs[[i]])
    all_distances[idx:(idx + length(distances) - 1)] <- distances
    idx <- idx + length(distances)
  }
  avg_d_G_Mi[i] <- mean(all_distances)
}

# Compute the sum part of the test statistic
sum_term <- 0
for (i in seq_along(populations)) {
  term1 <- (n_i[i] / (n_i[i] - 1)) * avg_d_Gi_Mi[i]
  term2 <- (n / (n - 1)) * avg_d_G_Mi[i]
  sum_term <- sum_term + sqrt(n_i[i]) * (term1 - term2)
}

a=1
# Compute the test statistic T
T_value <- (sqrt(m) / a) * sum_term
print(paste("Test statistic T_value:", T_value))


p_value_pnorm <- pnorm(T_value)
print(paste("Test statistic p_value_pnorm:", p_value_pnorm))


#--------------------#
#--------------------#
#--------------------#
#--------------------#
#--------------------#
#--------------------#
#--------------------#
#--------------------#



# Run function
result <- identify_critical_links(populations, alpha = 0.05, method = "fisher")

# Check results
expect_true(!is.null(result$critical_edges))

# Only one edge should be removed
expect_equal(nrow(result$critical_edges), 1)
expect_equal(length(result$edges_removed), 1)

# Check that edge (2,3) was identified and removed
edge_removed <- result$edges_removed[[1]]
edge_found <- ((edge_removed[1] == 2 && edge_removed[2] == 3) ||
                 (edge_removed[1] == 3 && edge_removed[2] == 2))
expect_true(edge_found)

# After removal, both populations should be identical
for (g in seq_along(result$modified_populations$A)) {
expect_equal(result$modified_populations$A[[g]], result$modified_populations$B[[g]])
}


























































#----------------------------------------------------------------#
# Test 20: OneEdgeOff_5x5
# Create base adjacency matrix with moderate density
set.seed(999)
base_mat <- matrix(rbinom(25, 1, 0.4), 5, 5)
diag(base_mat) <- 0
# Force symmetry
base_mat[lower.tri(base_mat)] <- base_mat[upper.tri(base_mat)]

# Population A: Force edge (1,4) = 1
A_mat <- base_mat
A_mat[1, 4] <- 1
A_mat[4, 1] <- 1  # Keep symmetry
A <- list(A_mat, A_mat)

# Population B: Force edge (1,4) = 0
B_mat <- base_mat
B_mat[1, 4] <- 0
B_mat[4, 1] <- 0  # Keep symmetry
B <- list(B_mat, B_mat)

# Create populations
populations <- list(A = A, B = B)

# Run function
result <- identify_critical_links(populations, alpha = 0.99, method = "fisher")

# Check results
expect_true(!is.null(result$critical_edges))

# Only one edge should be removed
expect_equal(nrow(result$critical_edges), 1)
expect_equal(length(result$edges_removed), 1)

# Check that edge (1,4) was identified and removed
edge_removed <- result$edges_removed[[1]]
edge_found <- ((edge_removed[1] == 1 && edge_removed[2] == 4) ||
                 (edge_removed[1] == 4 && edge_removed[2] == 1))
expect_true(edge_found)

# After removal, both populations should be identical
for (g in seq_along(result$modified_populations$A)) {
  expect_equal(result$modified_populations$A[[g]], result$modified_populations$B[[g]])
}

# The edge (1,4) should now be 0 in population A
expect_equal(result$modified_populations$A[[1]][1, 4], 0)
expect_equal(result$modified_populations$A[[1]][4, 1], 0)

#---------------- Another test ---------------#
# Test data for compute_test_statistic that should produce a highly negative T
test_negative_T_statistic <- function() {
  # Create two populations with fundamentally different structures

  # Population A: Very dense graphs
  pop_A <- list()
  for (i in 1:10) {
    # Create a dense graph with 90% of possible edges
    matrix_A <- matrix(1, nrow = 20, ncol = 20)
    diag(matrix_A) <- 0  # No self-loops
    # Add tiny random variations to avoid identical graphs
    random_indices <- sample(1:400, 10)
    matrix_A[random_indices] <- 0
    pop_A[[i]] <- matrix_A
  }

  # Population B: Very sparse graphs
  pop_B <- list()
  for (i in 1:10) {
    # Create a sparse graph with only 10% of possible edges
    matrix_B <- matrix(0, nrow = 20, ncol = 20)
    # Add tiny random variations
    random_indices <- sample(1:400, 40)
    matrix_B[random_indices] <- 1
    matrix_B <- matrix_B + t(matrix_B)  # Make symmetric
    matrix_B[matrix_B > 1] <- 1  # Ensure binary
    diag(matrix_B) <- 0  # No self-loops
    pop_B[[i]] <- matrix_B
  }

  # Combine populations
  populations <- list(A = pop_A, B = pop_B)

  # Compute test statistic
  T_value <- compute_test_statistic(populations, a = 1)
  return(T_value)
}

# Run the test
set.seed(123)  # For reproducibility
T_value <- test_negative_T_statistic()
print(paste("Test statistic T:", T_value))

if (m == 2) {
  # For 2 populations, be more sensitive
  critical_threshold <- 0.05
} else {
  # For 3+ populations, adjust threshold
  critical_threshold <- 0.07
}

# Calculate p-value (lower T = higher p-value)
p_value <- exp(-T_value / critical_threshold)

# Ensure p-value is between 0 and 1
p_value <- max(0, min(p_value, 1))
print(paste("Test statistic p_value:", p_value))


p_value_pnorm <- pnorm(T_value)
print(paste("Test statistic p_value_pnorm:", p_value_pnorm))



