#Install R packages for building annotation tables (if needed) and import them
ensure.packages <- function(cran.pkgs = character(0), bioc.pkgs = character(0)) {
  
  # CRAN
  missing.cran <- cran.pkgs[!sapply(cran.pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing.cran) > 0) {
    message("Installing missing CRAN packages: ", paste(missing.cran, collapse = ", "))
    install.packages(missing.cran, dependencies = TRUE, repos = "https://cloud.r-project.org")
  }
  
  # Bioconductor
  missing.bioc <- bioc.pkgs[!sapply(bioc.pkgs, requireNamespace, quietly = TRUE)]
  if (length(missing.bioc) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    message("Installing missing Bioconductor packages: ", paste(missing.bioc, collapse = ", "))
    BiocManager::install(missing.bioc, update = FALSE, ask = FALSE)
  }
  
  # Load packages
  invisible(sapply(c(cran.pkgs, bioc.pkgs), library, character.only = TRUE))
}
