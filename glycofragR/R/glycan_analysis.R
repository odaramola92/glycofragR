#' GlycanAnalysis wrapper
#'
#' High-level interface for glycan analysis combining structure prediction,
#' fragment generation, and Excel export.
#'
#' @param glycan_code 4-5 digit composition code
#' @param glycan_type Type of glycan: "N" or "O" (default: "N")
#' @param modification_type Reducing end modification (default: "permethylated_reduced")
#' @param max_structures Maximum number of structures to predict (default: 100)
#' @param isomer_sensitive Treat mirror images as distinct (default: FALSE)
#' @param fragment_types Vector of fragment types (default: c("BY"))
#' @param charges Vector of charge states (default: c(1, 2, 3))
#' @param visualize_structure Which structures to visualize: "all", "best", integer, or NULL
#' @param preferred_core For O-glycans, select specific core (1-8) or NULL
#' @return GlycanAnalysis object
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' analysis <- glycan_analysis_new(
#'   glycan_code = "4501",
#'   glycan_type = "N",
#'   modification_type = "permethylated_reduced"
#' )
#' }
glycan_analysis_new <- function(
  glycan_code,
  glycan_type = "N",
  modification_type = "permethylated_reduced",
  max_structures = 100,
  isomer_sensitive = FALSE,
  fragment_types = c("BY"),
  charges = c(1, 2, 3),
  visualize_structure = NULL,
  preferred_core = NULL
) {
  gf <- .glycofrag_py()
  gf$GlycanAnalysis(
    glycan_code = glycan_code,
    glycan_type = glycan_type,
    modification_type = modification_type,
    max_structures = as.integer(max_structures),
    isomer_sensitive = isomer_sensitive,
    fragment_types = as.list(fragment_types),
    charges = as.list(as.integer(charges)),
    visualize_structure = visualize_structure,
    preferred_core = preferred_core
  )
}

#' Export GlycanAnalysis to Excel
#'
#' @param analysis GlycanAnalysis object
#' @param output_file Path to output Excel file
#' @param include_prediction Include structure prediction sheet (default: TRUE)
#' @return NULL invisibly
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' analysis <- glycan_analysis_new("4501", "N")
#' glycan_analysis_export(analysis, "output.xlsx")
#' }
glycan_analysis_export <- function(analysis, output_file, include_prediction = TRUE) {
  analysis$export_to_excel(output_file, include_prediction = include_prediction)
  invisible(NULL)
}

#' Get GlycanAnalysis structures
#'
#' @param analysis GlycanAnalysis object
#' @return Data frame of predicted structures with structure_id, structure_str, and structure_repr
#' @details
#' Returns a data.frame containing information about all predicted structures.
#' Each row represents one predicted glycan structure.
#' @export
glycan_analysis_structures <- function(analysis) {
  # Get structures as R-compatible list format
  structures_list <- analysis$get_structures_as_list()
  # Convert to data.frame
  .list_dict_to_df(structures_list)
}

#' Get theoretical fragments from GlycanAnalysis
#'
#' @param analysis GlycanAnalysis object
#' @return Data frame of theoretical fragments
#' @details
#' Returns a native R data.frame with theoretical glycan fragments.
#' Each row represents one fragment ion variant.
#' Columns typically include: fragment_id, mass, composition, charge, m_z, etc.
#' @export
glycan_analysis_theoretical_df <- function(analysis) {
  # Get R-compatible list of dicts from Python
  frag_list <- analysis$get_theoretical_df_as_list()
  # Convert to data.frame
  .list_dict_to_df(frag_list)
}

#' Get clean theoretical fragments from GlycanAnalysis
#'
#' @param analysis GlycanAnalysis object
#' @return Data frame of unique/deduplicated fragments
#' @details
#' Returns a native R data.frame with deduplicated theoretical fragments.
#' Fragments are deduplicated by mass and composition.
#' This usually contains fewer rows than theoretical_df due to deduplication.
#' @export
glycan_analysis_clean_df <- function(analysis) {
  # Get R-compatible list of dicts from Python
  frag_list <- analysis$get_clean_df_as_list()
  # Convert to data.frame
  .list_dict_to_df(frag_list)
}

#' Get summary from GlycanAnalysis
#'
#' @param analysis GlycanAnalysis object
#' @return Data frame with analysis summary
#' @details
#' Returns a native R data.frame with summarized analysis information.
#' Includes glycan code, type, modification type, masses, and QC information.
#' @export
glycan_analysis_summary <- function(analysis) {
  # Get R-compatible list of dicts from Python
  summary_list <- analysis$get_summary_df_as_list()
  # Convert to data.frame
  .list_dict_to_df(summary_list)
}

#' Visualize structure from GlycanAnalysis
#'
#' @param analysis GlycanAnalysis object
#' @param structure_number Which structure to visualize (1-based)
#' @param output_path Path to save image
#' @param show Display the plot (default: FALSE)
#' @param show_node_numbers Show node numbers (default: TRUE)
#' @param ... Additional visualization parameters
#' @return NULL invisibly
#' @export
glycan_analysis_visualize <- function(
  analysis,
  structure_number,
  output_path = NULL,
  show = FALSE,
  show_node_numbers = TRUE,
  ...
) {
  args <- list(
    structure_number = as.integer(structure_number),
    output_path = output_path,
    show = show,
    show_node_numbers = show_node_numbers,
    ...
  )
  do.call(analysis$visualize_structure, args)
  invisible(NULL)
}
