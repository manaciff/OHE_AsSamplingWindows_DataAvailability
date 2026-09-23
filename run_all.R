# ==============================================================================
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
#
# Runs the whole chain in order.
#
# Each script verifies what it receives from the previous one and stops rather
# than analysing a reduced sample, so a failure here is informative: read the
# message before rerunning.
#
# Script 00 is the slow one, a few minutes. Script 05 runs 9,999 permutations
# twice and takes about ten. The rest are seconds to two minutes.
#
# Script 02 must run before script 03: it writes the Moran eigenvectors that the
# spatial sensitivity check of script 03 needs.
#
# WHERE THE SCRIPTS ARE LOOKED FOR
#   Earlier versions called source(here("R", s)). That fails whenever the active
#   RStudio project is the working folder rather than the repository, because
#   here() then resolves to the working folder, in which R/ does not exist:
#
#     cannot open file 'D:/Duda_Nacif_TCC/R/00_recompute_class_metrics.R'
#
#   The chain is now located from the position of this file, and here() is left
#   to do what it does well, which is to resolve Dados/ and Outputs/ inside the
#   working folder. If the location of this file cannot be determined, which
#   happens when the code is pasted into the console instead of sourced, a short
#   list of candidate folders is tried and the first one holding the chain is
#   used.
# ==============================================================================
library(here)

# Scripts 00 and 01 are commented out on purpose. Script 00 rebuilds the class
# metrics and script 01 fits the ordinations; both are memory intensive and
# script 01 ends the R session on the machine used here, so they are run one at
# a time from a clean session before this file. Uncomment them to reproduce the
# chain from the raw inputs. Everything below runs unattended.
#
# Script 08 is part of the chain, not an extra: it writes Table S42, the
# sampling effort bounds cited in the Results. It reads the model of Table 3 and
# therefore runs after script 05.
scripts <- c(
  #"00_recompute_class_metrics.R",
  #"01_prepare_data_and_ordination.R",
  "02_constrained_ordination_and_partitioning.R",
  "03_univariate_model_selection.R",
  "04_geometric_coupling_diagnostic.R",
  "05_ratio_artefact_test.R",
  "08_sampling_effort_sensitivity.R"
)

# ------------------------------------------------------------------------------
# Locate this file, then the folder holding the chain
# ------------------------------------------------------------------------------
this_file_path <- function() {
  # Rscript run_all.R
  args <- commandArgs(trailingOnly = FALSE)
  hit  <- grep("^--file=", args, value = TRUE)
  if (length(hit)) return(normalizePath(sub("^--file=", "", hit[1]), mustWork = FALSE))
  # source("run_all.R") from the console or from another script
  for (i in seq_len(sys.nframe())) {
    ofile <- sys.frame(i)$ofile
    if (!is.null(ofile)) return(normalizePath(ofile, mustWork = FALSE))
  }
  # the Source button of RStudio
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (!is.null(p) && nzchar(p)) return(normalizePath(p, mustWork = FALSE))
  }
  NA_character_
}

self <- this_file_path()

candidates <- c(
  if (!is.na(self)) file.path(dirname(self), "R"),
  file.path(getwd(), "R"),
  file.path(getwd(), "Repositorio_GitHub", "R"),
  here("R"),
  here("Repositorio_GitHub", "R")
)
candidates <- unique(candidates[!is.na(candidates)])

has_chain <- function(p) dir.exists(p) && all(file.exists(file.path(p, scripts)))
found     <- candidates[vapply(candidates, has_chain, logical(1))]

if (!length(found)) {
  stop("The six analysis scripts were not found. Folders tried:\n  ",
       paste(candidates, collapse = "\n  "),
       "\nOpen run_all.R in RStudio and press Source, or run it with ",
       "Rscript from the repository folder.", call. = FALSE)
}
script_dir <- found[1]

message("Scripts     : ", script_dir)
message("Project root: ", here())
message("Data        : ", here("Dados", "Processados"))
message("Outputs     : ", here("Outputs", "Manuscrito"))

# The anchor check exists because of a real failure. On 29 August 2026 the chain
# ran to completion and wrote every file into Repositorio_GitHub/Outputs instead of
# into the Outputs folder of the working directory. The cause is the .here file kept
# in this folder so that someone who clones only the repository can run the chain:
# it is an anchor for the here package, so whenever the R session starts inside
# Repositorio_GitHub, here() resolves to this folder and not to the folder that
# holds ANALISES_TCC.Rproj. Nothing errors, and the outputs land one level down.
if (!file.exists(file.path(here(), "ANALISES_TCC.Rproj"))) {
  stop("here() resolved the project root to\n  ", here(),
       "\nwhich does not contain ANALISES_TCC.Rproj, so the chain would write its ",
       "outputs there instead of into the working folder.\n",
       "Fix: open ANALISES_TCC.Rproj and press Source on this file, or set the ",
       "working directory to the folder that holds it before sourcing. The .here ",
       "file inside Repositorio_GitHub is what pulls the anchor down, and it is ",
       "kept so that the repository also runs on its own.", call. = FALSE)
}

for (p in c(here("Dados", "Processados"), here("Outputs", "Manuscrito"))) {
  if (!dir.exists(p)) {
    stop("Folder not found: ", p,
         "\nhere() resolved the project root to ", here(),
         ". Open ANALISES_TCC.Rproj, which is at the top of the working ",
         "folder, before running the chain.", call. = FALSE)
  }
}

# The log is written to disk and flushed line by line, so a session that dies
# mid-run still leaves a record of where it stopped.
.log_file <- file.path(here("Outputs", "Manuscrito"), "run_all_log.txt")
.say <- function(...) {
  txt <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  message(txt)
  cat(txt, "\n", file = .log_file, append = TRUE, sep = "")
}
cat("", file = .log_file)
.say("run_all started")
.say("R           : ", R.version.string)
.say("working dir : ", getwd())
.say("here()      : ", here())
.say("scripts     : ", script_dir)
.say("outputs     : ", here("Outputs", "Manuscrito"))

.stamp <- function() {
  f <- list.files(here("Outputs", "Manuscrito"), recursive = TRUE, full.names = TRUE)
  stats::setNames(file.mtime(f), f)
}

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------
.failed <- character(0)
for (.script in scripts) {
  message("\n", strrep("=", 78))
  .say("RUNNING  ", .script)
  message(strrep("=", 78))
  .before <- .stamp()
  .t0 <- Sys.time()
  ok <- tryCatch({
    source(file.path(script_dir, .script), echo = FALSE)
    TRUE
  }, error = function(e) {
    .say("ERROR in ", .script, ": ", conditionMessage(e))
    FALSE
  })
  .mins <- as.numeric(difftime(Sys.time(), .t0, units = "mins"))
  .after <- .stamp()
  .novos <- setdiff(names(.after), names(.before))
  .comuns <- intersect(names(.after), names(.before))
  .written <- length(.novos) + sum(.after[.comuns] > .before[.comuns], na.rm = TRUE)
  if (ok) {
    .say(sprintf("%s finished in %.1f min, %d files written", .script, .mins, .written))
    if (.written == 0) {
      .say("WARNING: ", .script, " wrote no file. Check here() and the output path.")
    }
  } else {
    .failed <- c(.failed, .script)
    .say(sprintf("%s stopped after %.1f min", .script, .mins))
  }
  .say(sprintf("memory in use: %.0f MB", sum(gc()[, 2])))
}
if (length(.failed)) {
  .say("chain incomplete. Scripts that stopped: ", paste(.failed, collapse = ", "))
} else {
  .say("all scripts finished")
}

# The closing note reads the sample size from the file the chain just wrote,
# rather than asserting it. A hand-maintained claim in the runner is one more
# thing that can fall out of step with the analysis.
raw <- file.path(here("Dados", "Processados"), "Data_Raw_WithCoords.csv")
n_units <- if (file.exists(raw)) nrow(read.csv(raw)) else NA_integer_
.say(sprintf("Chain complete. Sampling units in Data_Raw_WithCoords.csv: %s",
                ifelse(is.na(n_units), "file not found", n_units)))
message("The inferential estimates are in Table_S_PartIII_ConfirmatoryModel.csv")
message("and in TableS_PartV_Uncoupled_Model.csv. The dredge output describes")
message("selection uncertainty and is not a second set of estimates to test.")
