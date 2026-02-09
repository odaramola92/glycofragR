#' glycofragR
#'
#' R wrapper for the glycofrag Python package.
#'
#' @docType package
#' @name glycofragR
#' @importFrom reticulate import py_available py_to_r use_condaenv
NULL

# Declare global variables to avoid R CMD check NOTEs
utils::globalVariables(c(".glycofrag_py", ".glycofrag_visualizer"))
