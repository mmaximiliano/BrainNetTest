options(timeout = 600)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- sub("^--file=", "", script_arg)
  setwd(dirname(normalizePath(script_path)))
}

if (!requireNamespace("BrainNetTest", quietly = TRUE) ||
  as.character(utils::packageVersion("BrainNetTest")) != "1.0.0") {
  stop("Primero correr 00_install_brainnettest.R.")
}

pipeline_url <- paste0(
  "https://raw.githubusercontent.com/mmaximiliano/BrainNetTest/",
  "695f534/paper_source/data-application.R"
)
expected_md5 <- "0bbb3cd50e0323aaff51618c5fa0e254"

dir.create("pipeline", showWarnings = FALSE)
pipeline_file <- file.path(
  "pipeline",
  "data-application-695f534.R"
)

pipeline_is_valid <- file.exists(pipeline_file) &&
  identical(unname(tools::md5sum(pipeline_file)), expected_md5)
if (!pipeline_is_valid) {
  if (file.exists(pipeline_file)) {
    unlink(pipeline_file)
  }
  message("Descargando el pipeline clínico fijado al commit 695f534...")
  temporary_pipeline <- tempfile(
    pattern = "data-application-",
    tmpdir = "pipeline"
  )
  on.exit(unlink(temporary_pipeline), add = TRUE)
  utils::download.file(
    pipeline_url,
    temporary_pipeline,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )
  downloaded_md5 <- unname(tools::md5sum(temporary_pipeline))
  if (!identical(downloaded_md5, expected_md5)) {
    stop(
      "Checksum incorrecto para el pipeline descargado. Esperado: ",
      expected_md5, "; observado: ", downloaded_md5, "."
    )
  }
  if (!file.rename(temporary_pipeline, pipeline_file)) {
    stop("No se pudo promover el pipeline verificado.")
  }
}

observed_md5 <- unname(tools::md5sum(pipeline_file))
if (!identical(observed_md5, expected_md5)) {
  stop(
    "Checksum incorrecto para el pipeline. Esperado: ", expected_md5,
    "; observado: ", observed_md5,
    ". No se ejecutó código remoto."
  )
}

message("Checksum del pipeline verificado: ", observed_md5)

cache_dir <- Sys.getenv("ABIDE_CACHE", "abide_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
phenotypic_file <- file.path(
  cache_dir,
  "Phenotypic_V1_0b_preprocessed1.csv"
)
phenotypic_url <- paste0(
  "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/",
  "Phenotypic_V1_0b_preprocessed1.csv"
)
expected_phenotypic_md5 <- "033f8aac3da1066ff26fc52d027b3964"
phenotypic_is_valid <- file.exists(phenotypic_file) &&
  identical(
    unname(tools::md5sum(phenotypic_file)),
    expected_phenotypic_md5
  )
if (!phenotypic_is_valid) {
  if (file.exists(phenotypic_file)) {
    unlink(phenotypic_file)
  }
  temporary_phenotypic <- tempfile(
    pattern = "phenotypic-",
    tmpdir = cache_dir
  )
  on.exit(unlink(temporary_phenotypic), add = TRUE)
  utils::download.file(
    phenotypic_url,
    temporary_phenotypic,
    mode = "wb",
    method = "libcurl",
    quiet = FALSE
  )
  downloaded_md5 <- unname(tools::md5sum(temporary_phenotypic))
  if (!identical(downloaded_md5, expected_phenotypic_md5)) {
    stop(
      "Checksum incorrecto para el phenotypic file. Esperado: ",
      expected_phenotypic_md5, "; observado: ", downloaded_md5, "."
    )
  }
  if (!file.rename(temporary_phenotypic, phenotypic_file)) {
    stop("No se pudo promover el phenotypic file verificado.")
  }
}

existing_manifest_file <- "application/input_manifest.csv"
if (file.exists(existing_manifest_file)) {
  existing_manifest <- utils::read.csv(
    existing_manifest_file,
    stringsAsFactors = FALSE
  )
  if (!all(c("kind", "file", "bytes", "md5") %in%
    names(existing_manifest))) {
    stop("application/input_manifest.csv tiene un schema inválido.")
  }
  existing_paths <- vapply(seq_len(nrow(existing_manifest)), function(index) {
    switch(existing_manifest$kind[[index]],
      pipeline = pipeline_file,
      phenotypic = phenotypic_file,
      roi = file.path(
        cache_dir,
        "rois",
        existing_manifest$file[[index]]
      ),
      stop("Kind inválido en input_manifest.csv.")
    )
  }, character(1L))
  existing_files_present <- file.exists(existing_paths)
  existing_md5 <- rep(NA_character_, nrow(existing_manifest))
  existing_md5[existing_files_present] <- unname(tools::md5sum(
    existing_paths[existing_files_present]
  ))
  invalid <- !existing_files_present |
    existing_md5 != existing_manifest$md5
  if (any(invalid)) {
    stop(
      "El cache no coincide con application/input_manifest.csv: ",
      paste(existing_paths[invalid], collapse = ", "),
      ". Borrar los archivos inválidos y el manifest antes de regenerar."
    )
  }
}

phenotypic <- utils::read.csv(
  phenotypic_file,
  stringsAsFactors = FALSE
)
phenotypic$func_mean_fd <- suppressWarnings(
  as.numeric(phenotypic$func_mean_fd)
)
selected <- phenotypic[
  phenotypic$FILE_ID != "no_filename" &
    phenotypic$SITE_ID == "NYU" &
    !is.na(phenotypic$func_mean_fd) &
    phenotypic$func_mean_fd < 0.2,
]
if (nrow(selected) != 171L) {
  stop("El phenotypic filter no produjo las 171 personas esperadas.")
}

roi_dir <- file.path(cache_dir, "rois")
dir.create(roi_dir, recursive = TRUE, showWarnings = FALSE)
roi_base_url <- paste0(
  "https://s3.amazonaws.com/fcp-indi/data/Projects/ABIDE_Initiative/",
  "Outputs/cpac/filt_noglobal/rois_aal"
)
validate_roi <- function(path) {
  if (!file.exists(path) || file.size(path) == 0L) {
    return(FALSE)
  }
  parsed <- try(
    utils::read.table(path, header = TRUE, comment.char = ""),
    silent = TRUE
  )
  !inherits(parsed, "try-error") &&
    nrow(parsed) > 2L &&
    ncol(parsed) == 116L
}

message("Verificando/descargando ", nrow(selected), " ROI time series...")
for (file_id in selected$FILE_ID) {
  destination <- file.path(
    roi_dir,
    paste0(file_id, "_rois_aal.1D")
  )
  if (!validate_roi(destination)) {
    if (file.exists(destination)) {
      unlink(destination)
    }
    temporary_roi <- tempfile(
      pattern = paste0(file_id, "-"),
      tmpdir = roi_dir
    )
    url <- paste0(
      roi_base_url,
      "/",
      file_id,
      "_rois_aal.1D"
    )
    utils::download.file(
      url,
      temporary_roi,
      mode = "wb",
      method = "libcurl",
      quiet = TRUE
    )
    if (!validate_roi(temporary_roi)) {
      unlink(temporary_roi)
      stop("ROI download inválido: ", file_id)
    }
    if (!file.rename(temporary_roi, destination)) {
      unlink(temporary_roi)
      stop("No se pudo promover ROI verificado: ", file_id)
    }
  }
}

message("Ejecutando descarga y construcción de conectomas ABIDE...")
sys.source(pipeline_file, envir = globalenv())

required <- c(
  "application/abide_connectomes.rds",
  "application/abide_subjects.csv",
  "application/abide_provenance.csv"
)
missing <- required[!file.exists(required)]
if (length(missing) > 0L) {
  stop("Faltan outputs del pipeline: ", paste(missing, collapse = ", "))
}

subjects <- utils::read.csv(
  "application/abide_subjects.csv",
  stringsAsFactors = FALSE
)
provenance <- utils::read.csv(
  "application/abide_provenance.csv",
  stringsAsFactors = FALSE
)
connectomes <- readRDS("application/abide_connectomes.rds")
if (nrow(subjects) != 171L ||
  sum(subjects$group == "control") != 98L ||
  sum(subjects$group == "autism") != 73L ||
  sum(subjects$balanced & subjects$group == "control") != 55L ||
  sum(subjects$balanced & subjects$group == "autism") != 55L ||
  provenance$n_nodes != 114L ||
  provenance$n_edges != 6441L ||
  length(connectomes$node_labels) != 114L) {
  stop("Los outputs ABIDE no coinciden con el cohort/schema esperado.")
}

input_files <- c(
  pipeline_file,
  phenotypic_file,
  list.files(
    file.path(cache_dir, "rois"),
    pattern = "\\.1D$",
    full.names = TRUE
  )
)
input_manifest <- data.frame(
  kind = c(
    "pipeline",
    "phenotypic",
    rep("roi", length(input_files) - 2L)
  ),
  file = c(
    basename(pipeline_file),
    basename(phenotypic_file),
    basename(input_files[-c(1L, 2L)])
  ),
  bytes = file.info(input_files)$size,
  md5 = unname(tools::md5sum(input_files)),
  stringsAsFactors = FALSE
)
utils::write.csv(
  input_manifest,
  "application/input_manifest.csv",
  row.names = FALSE
)

message("Preparación ABIDE terminada. Outputs en application/.")
