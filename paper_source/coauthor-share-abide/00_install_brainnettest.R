options(repos = c(CRAN = "https://cloud.r-project.org"))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 1L) {
  script_path <- sub("^--file=", "", script_arg)
  setwd(dirname(normalizePath(script_path)))
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}

archive <- "BrainNetTest_1.0.0.tar.gz"
if (!file.exists(archive)) {
  stop("No se encontró ", archive, " en ", getwd(), ".")
}
expected_md5 <- "2df138ce0233276e45730a055ea59da9"
observed_md5 <- unname(tools::md5sum(archive))
if (!identical(observed_md5, expected_md5)) {
  stop(
    "El package archive no coincide con PACKAGE_ARCHIVE.txt. Esperado: ",
    expected_md5, "; observado: ", observed_md5, "."
  )
}

install.packages(archive, repos = NULL, type = "source")

installed <- as.character(utils::packageVersion("BrainNetTest"))
if (installed != "1.0.0") {
  stop("Se esperaba BrainNetTest 1.0.0; quedó instalada ", installed, ".")
}

message("BrainNetTest ", installed, " instalado correctamente.")
message("API principal: brainnet_data() -> brainnet_test() -> selected_edges().")
