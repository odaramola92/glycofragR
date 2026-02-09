#' Create a new Glycan object
#'
#' @param glycan_code 4-5 digit composition code (e.g., "4501" for HexNAc-4, Hex-5, Fuc-0, NeuAc-1)
#' @param glycan_type Type of glycan: "N" for N-glycan or "O" for O-glycan (default: "N")
#' @return Python Glycan object
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' glycan <- glycan_new("4501", "N")
#' }
glycan_new <- function(glycan_code, glycan_type = "N") {
  gf <- .glycofrag_py()
  gf$Glycan(glycan_code, glycan_type = glycan_type)
}

#' Predict glycan structures
#'
#' Predict possible structures from a glycan composition.
#'
#' @param glycan Glycan object created by glycan_new()
#' @return List of predicted structure graphs
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' glycan <- glycan_new("4501", "N")
#' structures <- glycan_predict_structures(glycan)
#' }
glycan_predict_structures <- function(glycan) {
  reticulate::py_to_r(glycan$predict_structures())
}

#' Generate glycan fragments
#'
#' Generate theoretical fragment ions for a glycan structure.
#'
#' @param glycan Glycan object
#' @param structure Structure graph from glycan_predict_structures()
#' @param modification_type Reducing end modification type (0-6, default: 0)
#' @param permethylated Logical; if TRUE, use permethylated masses (default: FALSE)
#' @param fragment_types Vector of fragment types (default: c("BY", "CZ"))
#' @return List containing fragments and metadata
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' glycan <- glycan_new("4501", "N")
#' structures <- glycan_predict_structures(glycan)
#' frags <- glycan_generate_fragments(glycan, structures[[1]])
#' }
glycan_generate_fragments <- function(
  glycan,
  structure,
  modification_type = 0,
  permethylated = FALSE,
  fragment_types = c("BY", "CZ")
) {
  reticulate::py_to_r(
    glycan$generate_fragments(
      structure,
      modification_type = modification_type,
      permethylated = permethylated,
      fragment_types = fragment_types
    )
  )
}

#' Calculate glycan mass
#'
#' Calculate mass from a glycan composition code.
#'
#' @param glycan_code 4-5 digit composition code
#' @return Numeric mass in Daltons
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' mass <- glycan_mass("4501")
#' }
glycan_mass <- function(glycan_code) {
  gf <- .glycofrag_py()
  calc <- gf$GlycanMassCalculator()
  calc$calculate_glycan_mass(glycan_code)
}

#' Visualize glycan structure
#'
#' Create a SNFG-compliant visualization of a glycan structure.
#'
#' @param structure Structure graph from glycan_predict_structures()
#' @param title Optional title for the plot
#' @param output_path Optional file path to save the image
#' @param show Logical; if TRUE, display the plot (default: FALSE)
#' @param ... Additional visualization parameters (vertical_gap, horizontal_spacing, node_size, figsize)
#' @return NULL invisibly
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' glycan <- glycan_new("4501", "N")
#' structures <- glycan_predict_structures(glycan)
#' glycan_visualize(structures[[1]], output_path = "glycan.png")
#' }
glycan_visualize <- function(structure, title = NULL, output_path = NULL, show = FALSE, ...) {
  viz <- .glycofrag_visualizer()$GlycanVisualizer
  args <- list(structure = structure, title = title, output_path = output_path, show = show, ...)
  do.call(viz$visualize, args)
}
