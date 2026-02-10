.onLoad <- function(libname, pkgname) {
  # Prevent reticulate from using ephemeral/auto environments
  if (Sys.getenv("GLYCOFRAG_DISABLE_AUTOENV") != "1") {
    # Try to use the glycofrag conda environment (non-blocking)
    tryCatch(
      reticulate::use_condaenv("glycofrag", required = FALSE),
      error = function(e) invisible(NULL)
    )
  }
  invisible(TRUE)
}
