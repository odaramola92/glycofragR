#' @keywords internal
.glycofrag_py <- function() {
  reticulate::import("glycofrag", delay_load = TRUE)
}

#' @keywords internal
.glycofrag_visualizer <- function() {
  reticulate::import("glycofrag.io.visualizer", delay_load = TRUE)
}

#' Set glycofrag conda environment
#'
#' Configure reticulate to use a specific conda environment containing glycofrag.
#'
#' @param env_name Name of the conda environment (default: "glycofrag")
#' @return TRUE invisibly
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' }
set_glycofrag_conda_env <- function(env_name = "glycofrag") {
  reticulate::use_condaenv(env_name, required = TRUE)
  invisible(TRUE)
}

#' List supported modifications
#'
#' Get a list of all supported peptide modifications.
#'
#' @param verbose Logical; if TRUE, print detailed information (default: TRUE)
#' @return List of supported modifications with masses and targets
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' mods <- list_supported_modifications()
#' }
list_supported_modifications <- function(verbose = TRUE) {
  gf <- .glycofrag_py()
  reticulate::py_to_r(gf$list_supported_modifications(verbose = verbose))
}
#' Convert Python list of dictionaries to R data.frame
#'
#' Internal helper function for converting Python list-of-dicts format
#' to R data.frame. Used by glycan_analysis and glycopeptide_analysis
#' functions to provide seamless R-compatible output.
#'
#' @param list_of_dicts Python list of dictionaries (from get_*_as_list() methods)
#' @return R data.frame, or empty data.frame if input is empty
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Internal use only
#' # Called by glycan_analysis_theoretical_df, glycan_analysis_clean_df, etc.
#' }
.list_dict_to_df <- function(list_of_dicts) {
  if (is.null(list_of_dicts) || length(list_of_dicts) == 0) {
    return(data.frame())
  }
  
  # Convert list of dicts to data.frame
  # Each element in the list becomes a row, preserving column names
  do.call(
    rbind,
    c(
      lapply(
        list_of_dicts,
        as.data.frame,
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      list(make.row.names = FALSE)
    )
  )
}