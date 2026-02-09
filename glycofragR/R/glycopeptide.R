#' Create a new Glycopeptide object
#'
#' @param peptide_sequence Amino acid sequence
#' @param glycan_code 4-5 digit glycan composition code
#' @param glycosylation_site Position of glycosylation (1-based indexing)
#' @param glycan_type Type of glycan: "N" or "O" (default: "N")
#' @param modifications Optional named list of peptide modifications with positions (1-based)
#' @param mod_string Optional modification string (e.g., "M:Ox; C:CAM")
#' @return Python Glycopeptide object
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' gp <- glycopeptide_new(
#'   peptide_sequence = "LCPDCPLLAPLNDSR",
#'   glycan_code = "4501",
#'   glycosylation_site = 12,
#'   glycan_type = "N"
#' )
#' }
glycopeptide_new <- function(
  peptide_sequence,
  glycan_code,
  glycosylation_site,
  glycan_type = "N",
  modifications = NULL,
  mod_string = NULL
) {
  gf <- .glycofrag_py()
  gf$Glycopeptide(
    peptide_sequence = peptide_sequence,
    glycan_code = glycan_code,
    glycosylation_site = glycosylation_site,
    glycan_type = glycan_type,
    modifications = modifications,
    mod_string = mod_string
  )
}

#' Generate glycopeptide fragments
#'
#' Generate all theoretical fragment ions for a glycopeptide.
#'
#' @param glycopeptide Glycopeptide object created by glycopeptide_new()
#' @return List of fragment ions (peptide, glycan, and glycopeptide fragments)
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' gp <- glycopeptide_new("LCPDCPLLAPLNDSR", "4501", 12, "N")
#' frags <- glycopeptide_fragments(gp)
#' }
glycopeptide_fragments <- function(glycopeptide) {
  reticulate::py_to_r(glycopeptide$generate_fragments())
}

#' Get glycopeptide summary
#'
#' Extract key information from a glycopeptide object.
#'
#' @param glycopeptide Glycopeptide object
#' @return List with peptide_sequence, glycan_code, glycosylation_site, glycan_type
#' @export
#' @examples
#' \dontrun{
#' set_glycofrag_conda_env("glycofrag")
#' gp <- glycopeptide_new("LCPDCPLLAPLNDSR", "4501", 12, "N")
#' summary <- glycopeptide_summary(gp)
#' }
glycopeptide_summary <- function(glycopeptide) {
  list(
    peptide_sequence = glycopeptide$peptide_sequence,
    glycan_code = glycopeptide$glycan_code,
    glycosylation_site = glycopeptide$glycosylation_site,
    glycan_type = glycopeptide$glycan_type
  )
}
