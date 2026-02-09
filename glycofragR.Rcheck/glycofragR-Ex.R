pkgname <- "glycofragR"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "glycofragR-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('glycofragR')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("glycan_analysis_export")
### * glycan_analysis_export

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glycan_analysis_export
### Title: Export GlycanAnalysis to Excel
### Aliases: glycan_analysis_export

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D analysis <- glycan_analysis_new("4501", "N")
##D glycan_analysis_export(analysis, "output.xlsx")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glycan_analysis_export", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("glycan_analysis_new")
### * glycan_analysis_new

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glycan_analysis_new
### Title: GlycanAnalysis wrapper
### Aliases: glycan_analysis_new

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D analysis <- glycan_analysis_new(
##D   glycan_code = "4501",
##D   glycan_type = "N",
##D   modification_type = "permethylated_reduced"
##D )
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glycan_analysis_new", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("glycan_generate_fragments")
### * glycan_generate_fragments

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glycan_generate_fragments
### Title: Generate glycan fragments
### Aliases: glycan_generate_fragments

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D glycan <- glycan_new("4501", "N")
##D structures <- glycan_predict_structures(glycan)
##D frags <- glycan_generate_fragments(glycan, structures[[1]])
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glycan_generate_fragments", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("glycan_mass")
### * glycan_mass

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glycan_mass
### Title: Calculate glycan mass
### Aliases: glycan_mass

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D mass <- glycan_mass("4501")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glycan_mass", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("glycan_new")
### * glycan_new

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glycan_new
### Title: Create a new Glycan object
### Aliases: glycan_new

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D glycan <- glycan_new("4501", "N")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glycan_new", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("glycan_predict_structures")
### * glycan_predict_structures

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glycan_predict_structures
### Title: Predict glycan structures
### Aliases: glycan_predict_structures

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D glycan <- glycan_new("4501", "N")
##D structures <- glycan_predict_structures(glycan)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glycan_predict_structures", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("glycan_visualize")
### * glycan_visualize

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glycan_visualize
### Title: Visualize glycan structure
### Aliases: glycan_visualize

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D glycan <- glycan_new("4501", "N")
##D structures <- glycan_predict_structures(glycan)
##D glycan_visualize(structures[[1]], output_path = "glycan.png")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glycan_visualize", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("glycopeptide_fragments")
### * glycopeptide_fragments

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glycopeptide_fragments
### Title: Generate glycopeptide fragments
### Aliases: glycopeptide_fragments

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D gp <- glycopeptide_new("LCPDCPLLAPLNDSR", "4501", 12, "N")
##D frags <- glycopeptide_fragments(gp)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glycopeptide_fragments", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("glycopeptide_new")
### * glycopeptide_new

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glycopeptide_new
### Title: Create a new Glycopeptide object
### Aliases: glycopeptide_new

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D gp <- glycopeptide_new(
##D   peptide_sequence = "LCPDCPLLAPLNDSR",
##D   glycan_code = "4501",
##D   glycosylation_site = 12,
##D   glycan_type = "N"
##D )
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glycopeptide_new", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("glycopeptide_summary")
### * glycopeptide_summary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: glycopeptide_summary
### Title: Get glycopeptide summary
### Aliases: glycopeptide_summary

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D gp <- glycopeptide_new("LCPDCPLLAPLNDSR", "4501", 12, "N")
##D summary <- glycopeptide_summary(gp)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("glycopeptide_summary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("list_supported_modifications")
### * list_supported_modifications

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: list_supported_modifications
### Title: List supported modifications
### Aliases: list_supported_modifications

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D mods <- list_supported_modifications()
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("list_supported_modifications", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("peptide_fragments")
### * peptide_fragments

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: peptide_fragments
### Title: Generate peptide fragments
### Aliases: peptide_fragments

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D peptide <- peptide_new("PEPTIDE")
##D frags <- peptide_fragments(peptide, charge_states = c(1, 2, 3))
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("peptide_fragments", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("peptide_new")
### * peptide_new

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: peptide_new
### Title: Create a new Peptide object
### Aliases: peptide_new

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
##D peptide <- peptide_new("PEPTIDE")
##D peptide_mod <- peptide_new("LCPDCPLLAPLNDSR", mod_string = "C:CAM")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("peptide_new", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("set_glycofrag_conda_env")
### * set_glycofrag_conda_env

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: set_glycofrag_conda_env
### Title: Set glycofrag conda environment
### Aliases: set_glycofrag_conda_env

### ** Examples

## Not run: 
##D set_glycofrag_conda_env("glycofrag")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("set_glycofrag_conda_env", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
