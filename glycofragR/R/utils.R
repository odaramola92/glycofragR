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

#' Install glycofrag Python Package
#'
#' One-time setup to install glycofrag and its dependencies into a 
#' dedicated conda environment. Run this once per system after installing 
#' glycofragR.
#'
#' @param envname Name of the conda environment to create (default: "glycofrag")
#' @param python_version Python version to use (default: "3.11")
#' @param force_reinstall Reinstall glycofrag even if already present (default: FALSE)
#'
#' @return NULL invisibly. Prints status messages.
#' @export
#' @details
#' This function:
#' 1. Installs or checks for Miniconda
#' 2. Creates a dedicated conda environment
#' 3. Installs glycofrag and dependencies
#'
#' After running this, use `library(glycofragR)` normally; the package
#' will automatically find the glycofrag environment.
#'
#' @examples
#' \dontrun{
#' glycofrag_install()
#' }
glycofrag_install <- function(
  envname = "glycofrag",
  python_version = "3.11",
  force_reinstall = FALSE
) {
  # Step 1: Ensure Miniconda is available
  envs <- tryCatch(reticulate::conda_list(), error = function(e) NULL)
  if (is.null(envs)) {
    message("Installing Miniconda...")
    reticulate::install_miniconda()
    envs <- reticulate::conda_list()
  }
  
  # Step 2: Check if environment already exists
  env_exists <- envname %in% envs$name
  
  # Step 3: Create environment if needed
  if (!env_exists) {
    message("Creating conda environment: ", envname)
    reticulate::conda_create(
      envname,
      packages = paste0("python=", python_version),
      channel = "conda-forge"
    )
  } else if (!force_reinstall) {
    message("✓ Conda environment '", envname, "' already exists")
  }
  
  # Step 4: Install core dependencies via conda
  message("Installing core dependencies...")
  reticulate::conda_install(
    envname,
    packages = c("pip", "numpy", "pandas", "matplotlib", "scipy", "scikit-learn"),
    channel = "conda-forge"
  )
  
  # Step 5: Install glycofrag via pip from TestPyPI
  message("Installing glycofrag from TestPyPI...")
  result <- system2(
    reticulate::conda_binary(),
    args = c(
      "run", "-n", envname, "pip", "install",
      "--index-url", "https://test.pypi.org/simple/",
      "--extra-index-url", "https://pypi.org/simple/",
      if (force_reinstall) "--force-reinstall" else NULL,
      "glycofrag"
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  
  if (!is.null(attr(result, "status")) && attr(result, "status") != 0) {
    stop("Failed to install glycofrag from TestPyPI", call. = FALSE)
  }
  
  # Step 6: Verify installation
  message("Verifying installation...")
  tryCatch(
    {
      reticulate::use_condaenv(envname, required = TRUE)
      gf <- reticulate::import("glycofrag", delay_load = FALSE)
      # Check that module was loaded successfully
      if (!is.null(gf)) {
        message("\n✓ Success! glycofrag installed in environment: ", envname)
        message("You can now use: library(glycofragR)")
      }
    },
    error = function(e) {
      stop(
        "✗ Installation verification failed.\n",
        "Error: ", conditionMessage(e), "\n",
        "Try running: glycofrag_install(force_reinstall = TRUE)",
        call. = FALSE
      )
    }
  )
  
  invisible(NULL)
}
