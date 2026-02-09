#' Create a new Peptide object
#'
#' @param sequence Amino acid sequence (single-letter code)
#' @param modifications Optional named list of modifications with positions (1-based)
#' @param mod_string Optional modification string (e.g., "M:Ox; S3:Phos")
#' @return Python Peptide object
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' peptide <- peptide_new("PEPTIDE")
#' peptide_mod <- peptide_new("LCPDCPLLAPLNDSR", mod_string = "C:CAM")
#' }
peptide_new <- function(sequence, modifications = NULL, mod_string = NULL) {
  gf <- .glycofrag_py()
  gf$Peptide(sequence, modifications = modifications, mod_string = mod_string)
}

#' Generate peptide fragments
#'
#' Generate theoretical peptide fragment ions (b, y, c, z).
#'
#' @param peptide Peptide object created by peptide_new()
#' @param charge_states Vector of charge states to generate (default: c(1, 2))
#' @return List of fragment ions
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' peptide <- peptide_new("PEPTIDE")
#' frags <- peptide_fragments(peptide, charge_states = c(1, 2, 3))
#' }
peptide_fragments <- function(peptide, charge_states = c(1, 2)) {
  reticulate::py_to_r(peptide$generate_all_fragments(charge_states = charge_states))
}
