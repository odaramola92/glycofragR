# Install required packages for GlycofragR wrapper

message("Installing required R packages...")

packages_to_install <- c("reticulate", "writexl", "dplyr")

for (pkg in packages_to_install) {
  if (!require(pkg, character.only = TRUE)) {
    message(paste("Installing", pkg, "..."))
    install.packages(pkg, repos = "https://cran.r-project.org")
  } else {
    message(paste(pkg, "already installed."))
  }
}

message("\n✓ Package installation complete!")
message("You can now run: source('GlycofragR_test_UPDATED.R')")
