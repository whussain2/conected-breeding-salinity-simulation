################################################################################
## 04_run_strategies.R -- additional established strategies simulated under the
## SAME parameterisation, so that the comparative entries in Supplemental
## Table S5 are outputs of this simulation rather than unattributed estimates.
##
##   PSRS  : closed-pool recurrent phenotypic selection
##   MAGIC : one-time broad-base multi-parent population, then closed genomic
##           recurrent selection (no recurrent donor inflow)
##
## Usage: Rscript R/04_run_strategies.R <nRep> <firstSeed> <outfile>
################################################################################

Sys.setenv(OMP_NUM_THREADS = "1")

## ------------------------------------------------------------------ paths --
## Locate the bundle this script belongs to, so it works whether it is run with
## Rscript, sourced, or run from the RStudio editor, and whatever the working
## directory happens to be.  The script's OWN folder is tried first, so if two
## copies of the bundle exist (for example a newer one unzipped inside an older
## one) each copy uses its own data.  Set CB_SIM_ROOT to override.
.cb_script_dir <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f  <- grep("^--file=", ca, value = TRUE)
  if (length(f)) return(dirname(sub("^--file=", "", f[1])))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- tryCatch(get0("ofile", envir = sys.frame(i), inherits = FALSE),
                   error = function(e) NULL)
    if (is.character(of) && length(of) == 1L) return(dirname(of))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path,
                  error = function(e) "")
    if (is.character(p) && nzchar(p)) return(dirname(p))
  }
  NA_character_
}

.cb_is_root <- function(p) file.exists(file.path(p, "R", "01_functions.R"))

.cb_candidates <- function() {
  cand <- character(0)
  env <- Sys.getenv("CB_SIM_ROOT")
  if (nzchar(env)) cand <- c(cand, env)
  sd <- .cb_script_dir()
  if (!is.na(sd)) cand <- c(cand, file.path(sd, ".."), sd)
  cand <- c(cand, getwd(), file.path(getwd(), ".."),
            file.path(getwd(), "..", ".."))
  cand <- unique(normalizePath(cand, mustWork = FALSE))
  cand[vapply(cand, .cb_is_root, logical(1))]
}

.cb_roots <- .cb_candidates()
if (!length(.cb_roots)) {
  stop("Could not locate the simulation bundle.\n",
       "  Expected a folder containing R/01_functions.R at or above: ", getwd(),
       "\n  Run this script from the bundle folder, e.g.\n",
       "      setwd('/path/to/connected-breeding-salinity-simulation')\n",
       "      source('R/08_figures.R')\n",
       "  or set it explicitly:\n",
       "      Sys.setenv(CB_SIM_ROOT = '/path/to/connected-breeding-salinity-simulation')",
       call. = FALSE)
}
ROOT   <- .cb_roots[1]
OUTDIR <- file.path(ROOT, "out")
FIGDIR <- file.path(ROOT, "figures")
RDIR   <- file.path(ROOT, "R")
message("bundle root: ", ROOT)
## ---------------------------------------------------------------------------

suppressMessages(source(file.path(RDIR, "01_functions.R")))

args    <- commandArgs(trailingOnly = TRUE)
nRep    <- if (length(args) >= 1) as.integer(args[1]) else 10
seed0   <- if (length(args) >= 2) as.integer(args[2]) else 9000
outfile <- if (length(args) >= 3) args[3] else file.path(OUTDIR, "raw_strategies.csv")

if (!grepl("^([A-Za-z]:)?[/\\\\]", outfile))
  outfile <- file.path(ROOT, outfile)
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

base <- defaultParams()
STRATS <- c("PSRS", "MAGIC")

first <- TRUE
for (r in seq_len(nRep)) {
  seed <- seed0 + r
  f <- makeFounders(base, seed = seed)
  donors <- drawDonors(f, base, 0)
  for (st in STRATS) {
    res <- try(runProgram(st, f, base, donors = donors, seed = seed),
               silent = TRUE)
    if (inherits(res, "try-error")) {
      message("FAILED ", st, " rep ", r, ": ", as.character(res)); next
    }
    res$scenario <- st
    res$axis     <- "strategy_comparison"
    res$level    <- st
    res$rep      <- r
    res$seed     <- seed
    write.table(res, outfile, sep = ",", row.names = FALSE,
                col.names = first, append = !first)
    first <- FALSE
  }
  message(sprintf("rep %d/%d done", r, nRep))
}
message("finished -> ", outfile)
