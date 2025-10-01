#!/usr/bin/env Rscript

options(repos = c(
  CRAN  = "https://cloud.r-project.org",                  # a globally‑replicated CRAN mirror
  BioC  = "https://bioconductor.org/packages/3.18/bioc"    # (optional) Bioconductor mirror
))

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# List of packages to install
packages_to_install <- c(
    "MEIGOR",
    "CellNOptR",
    "BoolNet",
    "here",
    "optparse",
    "dplyr",
    "tidyr",
    "readr",
    "data.table"
)

cat("Installing packages...\n")
for (pkg in packages_to_install) {
    cat(sprintf("Installing %s...\n", pkg))
    tryCatch({
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
        cat(sprintf("✓ %s installed successfully\n", pkg))
    }, error = function(e) {
        cat(sprintf("✗ Failed to install %s: %s\n", pkg, e$message))
    })
}

cat("Testing CellNOptR installation...\n")
tryCatch({
    library(CellNOptR)
    cat("✓ CellNOptR loaded successfully\n")
    cat(sprintf("CellNOptR version: %s\n", packageVersion("CellNOptR")))
}, error = function(e) {
    cat(sprintf("✗ Failed to load CellNOptR: %s\n", e$message))
})

cat("Installation complete!\n")
