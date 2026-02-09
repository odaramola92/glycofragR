.require_py <- function() {
  if (!reticulate::py_available(initialize = FALSE)) {
    stop("Python is not available. Use set_glycofrag_conda_env() first.", call. = FALSE)
  }
  invisible(TRUE)
}
