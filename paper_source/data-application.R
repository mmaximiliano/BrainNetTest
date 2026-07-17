###########################################################################
## Real-data application for:
## BrainNetTest: Global and Edge-Wise Inference for Binary Brain Network
## Populations
##
## This standalone script downloads openly shared, preprocessed
## resting-state fMRI derivatives from the Autism Brain Imaging Data
## Exchange (ABIDE I) and turns them into a population of aligned binary,
## undirected brain networks: one network per participant, in two
## independent groups (autism and control). It saves only the derived
## binary connectomes and provenance under `application/`. The manuscript
## reads those derived objects; it never downloads data at compile time.
##
## Data source and attribution (please cite if you reuse the derivatives):
##   * ABIDE I: Di Martino et al. (2014), doi:10.1038/mp.2013.78.
##   * Preprocessed derivatives (Neuro Bureau / Preprocessed Connectomes
##     Project): Craddock et al. (2013).
##   * AAL atlas: Tzourio-Mazoyer et al. (2002), doi:10.1006/nimg.2001.0978.
## ABIDE derivatives are shared for non-commercial research use under a
## Creative Commons Attribution-NonCommercial-ShareAlike license. The
## derived binary connectomes saved here inherit that license.
##
## Requirements: internet access and about 35 MB of downloads on the first
## run. Downloads are cached in `abide_cache/` and reused afterwards.
###########################################################################

options(stringsAsFactors = FALSE, timeout = 600)
suppressPackageStartupMessages(library("BrainNetTest"))

## ---- Analysis choices ----
site <- "NYU"                 # single site avoids multi-site confounds
pipeline <- "cpac"            # Configurable Pipeline for the Analysis of Connectomes
strategy <- "filt_noglobal"   # band-pass filtered, no global signal regression
derivative <- "rois_aal"      # AAL atlas mean time series (116 regions)
max_mean_fd <- 0.2            # motion quality control on mean framewise displacement
density <- 0.10               # proportional edge-density threshold per subject
age_bin_width <- 5             # years, for the balanced sensitivity sample
fd_bin_width <- 0.05           # mean FD, for the balanced sensitivity sample

cache_dir <- Sys.getenv("ABIDE_CACHE", "abide_cache")
roi_dir <- file.path(cache_dir, "rois")
output_dir <- "application"
dir.create(roi_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, showWarnings = FALSE)

## AAL region codes (as they appear in the .1D column headers) mapped to the
## anatomical names of Tzourio-Mazoyer et al. (2002), used for readable labels.
aal_names <- c(
  "2001" = "Precentral_L", "2002" = "Precentral_R",
  "2101" = "Frontal_Sup_L", "2102" = "Frontal_Sup_R",
  "2111" = "Frontal_Sup_Orb_L", "2112" = "Frontal_Sup_Orb_R",
  "2201" = "Frontal_Mid_L", "2202" = "Frontal_Mid_R",
  "2211" = "Frontal_Mid_Orb_L", "2212" = "Frontal_Mid_Orb_R",
  "2301" = "Frontal_Inf_Oper_L", "2302" = "Frontal_Inf_Oper_R",
  "2311" = "Frontal_Inf_Tri_L", "2312" = "Frontal_Inf_Tri_R",
  "2321" = "Frontal_Inf_Orb_L", "2322" = "Frontal_Inf_Orb_R",
  "2331" = "Rolandic_Oper_L", "2332" = "Rolandic_Oper_R",
  "2401" = "Supp_Motor_Area_L", "2402" = "Supp_Motor_Area_R",
  "2501" = "Olfactory_L", "2502" = "Olfactory_R",
  "2601" = "Frontal_Sup_Medial_L", "2602" = "Frontal_Sup_Medial_R",
  "2611" = "Frontal_Med_Orb_L", "2612" = "Frontal_Med_Orb_R",
  "2701" = "Rectus_L", "2702" = "Rectus_R",
  "3001" = "Insula_L", "3002" = "Insula_R",
  "4001" = "Cingulum_Ant_L", "4002" = "Cingulum_Ant_R",
  "4011" = "Cingulum_Mid_L", "4012" = "Cingulum_Mid_R",
  "4021" = "Cingulum_Post_L", "4022" = "Cingulum_Post_R",
  "4101" = "Hippocampus_L", "4102" = "Hippocampus_R",
  "4111" = "ParaHippocampal_L", "4112" = "ParaHippocampal_R",
  "4201" = "Amygdala_L", "4202" = "Amygdala_R",
  "5001" = "Calcarine_L", "5002" = "Calcarine_R",
  "5011" = "Cuneus_L", "5012" = "Cuneus_R",
  "5021" = "Lingual_L", "5022" = "Lingual_R",
  "5101" = "Occipital_Sup_L", "5102" = "Occipital_Sup_R",
  "5201" = "Occipital_Mid_L", "5202" = "Occipital_Mid_R",
  "5301" = "Occipital_Inf_L", "5302" = "Occipital_Inf_R",
  "5401" = "Fusiform_L", "5402" = "Fusiform_R",
  "6001" = "Postcentral_L", "6002" = "Postcentral_R",
  "6101" = "Parietal_Sup_L", "6102" = "Parietal_Sup_R",
  "6201" = "Parietal_Inf_L", "6202" = "Parietal_Inf_R",
  "6211" = "SupraMarginal_L", "6212" = "SupraMarginal_R",
  "6221" = "Angular_L", "6222" = "Angular_R",
  "6301" = "Precuneus_L", "6302" = "Precuneus_R",
  "6401" = "Paracentral_Lobule_L", "6402" = "Paracentral_Lobule_R",
  "7001" = "Caudate_L", "7002" = "Caudate_R",
  "7011" = "Putamen_L", "7012" = "Putamen_R",
  "7021" = "Pallidum_L", "7022" = "Pallidum_R",
  "7101" = "Thalamus_L", "7102" = "Thalamus_R",
  "8101" = "Heschl_L", "8102" = "Heschl_R",
  "8111" = "Temporal_Sup_L", "8112" = "Temporal_Sup_R",
  "8121" = "Temporal_Pole_Sup_L", "8122" = "Temporal_Pole_Sup_R",
  "8201" = "Temporal_Mid_L", "8202" = "Temporal_Mid_R",
  "8211" = "Temporal_Pole_Mid_L", "8212" = "Temporal_Pole_Mid_R",
  "8301" = "Temporal_Inf_L", "8302" = "Temporal_Inf_R",
  "9001" = "Cerebelum_Crus1_L", "9002" = "Cerebelum_Crus1_R",
  "9011" = "Cerebelum_Crus2_L", "9012" = "Cerebelum_Crus2_R",
  "9021" = "Cerebelum_3_L", "9022" = "Cerebelum_3_R",
  "9031" = "Cerebelum_4_5_L", "9032" = "Cerebelum_4_5_R",
  "9041" = "Cerebelum_6_L", "9042" = "Cerebelum_6_R",
  "9051" = "Cerebelum_7b_L", "9052" = "Cerebelum_7b_R",
  "9061" = "Cerebelum_8_L", "9062" = "Cerebelum_8_R",
  "9071" = "Cerebelum_9_L", "9072" = "Cerebelum_9_R",
  "9081" = "Cerebelum_10_L", "9082" = "Cerebelum_10_R",
  "9100" = "Vermis_1_2", "9110" = "Vermis_3", "9120" = "Vermis_4_5",
  "9130" = "Vermis_6", "9140" = "Vermis_7", "9150" = "Vermis_8",
  "9160" = "Vermis_9", "9170" = "Vermis_10"
)

base_url <- paste0(
  "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/Outputs/",
  pipeline, "/", strategy, "/", derivative
)
phenotypic_url <- paste0(
  "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/",
  "Phenotypic_V1_0b_preprocessed1.csv"
)

## ---- Select a fixed, reproducible subject subset ----
phenotypic_file <- file.path(cache_dir, "Phenotypic_V1_0b_preprocessed1.csv")
if (!file.exists(phenotypic_file)) {
  utils::download.file(phenotypic_url, phenotypic_file, quiet = TRUE)
}
phenotypic <- utils::read.csv(phenotypic_file)
phenotypic <- phenotypic[phenotypic$FILE_ID != "no_filename", ]
phenotypic$func_mean_fd <- suppressWarnings(as.numeric(phenotypic$func_mean_fd))
phenotypic <- phenotypic[
  phenotypic$SITE_ID == site &
    !is.na(phenotypic$func_mean_fd) &
    phenotypic$func_mean_fd < max_mean_fd,
]
## Deterministic order: group, then anonymized file identifier.
phenotypic <- phenotypic[order(phenotypic$DX_GROUP, phenotypic$FILE_ID), ]
phenotypic$group <- ifelse(phenotypic$DX_GROUP == 1L, "autism", "control")
phenotypic$graph_index <- ave(
  seq_len(nrow(phenotypic)),
  phenotypic$group,
  FUN = seq_along
)

## The full cohort has different sex and motion distributions between groups.
## Build a deterministic coarsened-balance sensitivity sample without using
## any network outcome: exact sex strata, five-year age bins, and 0.05 mean-FD
## bins. Within each stratum, anonymized file-ID order breaks ties.
phenotypic$balance_stratum <- interaction(
  phenotypic$SEX,
  floor(phenotypic$AGE_AT_SCAN / age_bin_width),
  floor(phenotypic$func_mean_fd / fd_bin_width),
  drop = TRUE
)
phenotypic$balanced <- FALSE
phenotypic$balance_id <- NA_integer_
next_balance_id <- 0L
for (stratum in sort(unique(as.character(phenotypic$balance_stratum)))) {
  autism_rows <- which(
    phenotypic$balance_stratum == stratum &
      phenotypic$group == "autism"
  )
  control_rows <- which(
    phenotypic$balance_stratum == stratum &
      phenotypic$group == "control"
  )
  n_pairs <- min(length(autism_rows), length(control_rows))
  if (n_pairs > 0L) {
    autism_rows <- autism_rows[seq_len(n_pairs)]
    control_rows <- control_rows[seq_len(n_pairs)]
    ids <- next_balance_id + seq_len(n_pairs)
    phenotypic$balanced[c(autism_rows, control_rows)] <- TRUE
    phenotypic$balance_id[autism_rows] <- ids
    phenotypic$balance_id[control_rows] <- ids
    next_balance_id <- next_balance_id + n_pairs
  }
}
if (next_balance_id < 2L) {
  stop("The coarsened-balance sensitivity sample has fewer than 2 pairs.")
}
message(
  "Selected ", nrow(phenotypic), " ", site, " participants (",
  sum(phenotypic$group == "control"), " control, ",
  sum(phenotypic$group == "autism"), " autism); balanced sensitivity sample: ",
  next_balance_id, " per group."
)

## ---- Download the per-subject region time series (cached) ----
for (file_id in phenotypic$FILE_ID) {
  destination <- file.path(roi_dir, paste0(file_id, "_", derivative, ".1D"))
  if (!file.exists(destination) || file.size(destination) == 0L) {
    url <- paste0(base_url, "/", file_id, "_", derivative, ".1D")
    utils::download.file(url, destination, quiet = TRUE)
  }
}

read_time_series <- function(file_id) {
  path <- file.path(roi_dir, paste0(file_id, "_", derivative, ".1D"))
  as.matrix(utils::read.table(path, header = TRUE, comment.char = ""))
}

## ---- Drop regions with no signal in any subject, keep a common node set ----
region_codes <- sub("^X\\.", "", colnames(read_time_series(phenotypic$FILE_ID[1L])))
zero_variance <- logical(length(region_codes))
for (file_id in phenotypic$FILE_ID) {
  zero_variance <- zero_variance |
    (apply(read_time_series(file_id), 2L, stats::sd) == 0)
}
kept <- which(!zero_variance)
node_labels <- unname(aal_names[region_codes[kept]])
if (anyNA(node_labels)) {
  stop("Unmapped AAL region code(s): ",
    paste(region_codes[kept][is.na(node_labels)], collapse = ", "))
}

## ---- Build one binary, undirected connectome per subject ----
## Pearson correlation, then a proportional threshold that keeps the
## strongest `density` fraction of edges, so every subject network has the
## same edge count and the binary contract is satisfied.
build_binary_connectome <- function(file_id) {
  series <- read_time_series(file_id)[, kept, drop = FALSE]
  correlation <- stats::cor(series)
  diag(correlation) <- 0
  upper <- upper.tri(correlation)
  threshold <- stats::quantile(correlation[upper], probs = 1 - density, type = 7)
  adjacency <- matrix(0L, ncol(correlation), ncol(correlation))
  upper_edges <- integer(sum(upper))
  upper_edges[correlation[upper] > threshold] <- 1L
  adjacency[upper] <- upper_edges
  adjacency <- adjacency + t(adjacency)
  diag(adjacency) <- 0L
  storage.mode(adjacency) <- "integer"
  adjacency
}

populations <- list(
  control = lapply(
    phenotypic$FILE_ID[phenotypic$group == "control"],
    build_binary_connectome
  ),
  autism = lapply(
    phenotypic$FILE_ID[phenotypic$group == "autism"],
    build_binary_connectome
  )
)
balanced_populations <- lapply(names(populations), function(group) {
  indices <- phenotypic$graph_index[
    phenotypic$group == group & phenotypic$balanced
  ]
  populations[[group]][indices]
})
names(balanced_populations) <- names(populations)

## Validate against the package contract before saving.
networks <- brainnet_data(populations, node_labels = node_labels)
balanced_networks <- brainnet_data(
  balanced_populations,
  node_labels = node_labels
)

connectomes <- list(
  populations = populations,
  balanced_populations = balanced_populations,
  node_labels = node_labels,
  meta = list(
    dataset = "ABIDE I (Preprocessed Connectomes Project)",
    site = site,
    pipeline = pipeline,
    strategy = strategy,
    derivative = derivative,
    atlas = "AAL",
    n_nodes = length(node_labels),
    dropped_regions = region_codes[zero_variance],
    max_mean_fd = max_mean_fd,
    density = density,
    balance = list(
      method = "coarsened 1:1 sensitivity sample",
      exact = "sex",
      age_bin_width = age_bin_width,
      fd_bin_width = fd_bin_width,
      n_per_group = next_balance_id
    ),
    license = "CC BY-NC-SA (ABIDE / Neuro Bureau derivatives)"
  )
)
saveRDS(connectomes, file.path(output_dir, "abide_connectomes.rds"), compress = "xz")

utils::write.csv(
  transform(
    phenotypic[, c(
      "FILE_ID", "group", "AGE_AT_SCAN", "SEX", "func_mean_fd",
      "balanced", "balance_id", "balance_stratum"
    )],
    balance_stratum = as.character(balance_stratum)
  ),
  file.path(output_dir, "abide_subjects.csv"),
  row.names = FALSE
)

provenance <- data.frame(
  dataset = "ABIDE I (Preprocessed Connectomes Project)",
  site = site,
  pipeline = pipeline,
  strategy = strategy,
  derivative = derivative,
  atlas = "AAL",
  n_control = sum(phenotypic$group == "control"),
  n_autism = sum(phenotypic$group == "autism"),
  n_balanced_control = sum(
    phenotypic$group == "control" & phenotypic$balanced
  ),
  n_balanced_autism = sum(
    phenotypic$group == "autism" & phenotypic$balanced
  ),
  n_nodes = length(node_labels),
  n_edges = length(node_labels) * (length(node_labels) - 1L) / 2L,
  max_mean_fd = max_mean_fd,
  density = density,
  age_bin_width = age_bin_width,
  fd_bin_width = fd_bin_width,
  package_version = as.character(utils::packageVersion("BrainNetTest")),
  r_version = R.version.string,
  connectomes_md5 = unname(tools::md5sum(
    file.path(output_dir, "abide_connectomes.rds")
  )),
  downloaded = as.character(Sys.Date())
)
utils::write.csv(
  provenance,
  file.path(output_dir, "abide_provenance.csv"),
  row.names = FALSE
)

message(
  "Saved ", provenance$n_control + provenance$n_autism,
  " binary connectomes (", provenance$n_nodes, " nodes, ",
  provenance$n_edges, " possible edges), including a balanced ",
  provenance$n_balanced_control, " + ", provenance$n_balanced_autism,
  " sensitivity sample, to ",
  file.path(output_dir, "abide_connectomes.rds"), "."
)
