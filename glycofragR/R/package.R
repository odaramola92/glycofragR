#' glycofragR
#'
#' R wrapper for the glycofrag Python package.
#'
#' @docType package
#' @name glycofragR
#' @importFrom reticulate import py_available py_to_r use_condaenv conda_list 
#'   conda_create conda_install conda_binary_exists py_install
NULL

# Declare global variables to avoid R CMD check NOTEs
utils::globalVariables(c(".glycofrag_py", ".glycofrag_visualizer"))
