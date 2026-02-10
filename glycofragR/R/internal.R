.require_py <- function() {
  if (!reticulate::py_available(initialize = FALSE)) {
    stop(
      call. = FALSE,
      "\n✗ Python/glycofrag is not available.\n\n",
      "To set up glycofrag one time, run:\n",
      "  glycofrag_install()\n\n",
      "Alternatively, configure an existing conda env:\n",
      "  set_glycofrag_conda_env(\"your_env_name\")"
    )
  }
  
  # Test if glycofrag module is actually available
  tryCatch(
    reticulate::import("glycofrag", delay_load = FALSE),
    error = function(e) {
      stop(
        call. = FALSE,
        "\n✗ Python package 'glycofrag' is not installed.\n\n",
        "To set it up, run:\n",
        "  glycofrag_install()\n\n",
        "Details: ", conditionMessage(e)
      )
    }
  )
  
  invisible(TRUE)
}
