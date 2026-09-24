################################################################################
## 03_run.R -- driver: runs every scenario across replicates and writes
##             out/raw_results.csv
##
## Usage:  Rscript R/03_run.R <nRep> <firstSeed> <outfile> [group]
##         group is optional: "core" or "peripheral"
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

suppressMessages({
  source(file.path(RDIR, "01_functions.R"))
  source(file.path(RDIR, "02_scenarios.R"))
})

args    <- commandArgs(trailingOnly = TRUE)
nRep    <- if (length(args) >= 1) as.integer(args[1]) else 20
seed0   <- if (length(args) >= 2) as.integer(args[2]) else 1000
outfile <- if (length(args) >= 3) args[3] else file.path(OUTDIR, "raw_results.csv")
group   <- if (length(args) >= 4) args[4] else NULL

## normalise the output path: a bare or relative name goes into the bundle's
## out/ folder, an absolute path is used as given
if (!grepl("^([A-Za-z]:)?[/\\\\]", outfile))
  outfile <- file.path(ROOT, outfile)
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

S    <- buildScenarios(group)
base <- defaultParams()

## G x E levels that need their own founder object (trait must be created with
## the appropriate varGxE at SimParam time)
gxeLevels <- sort(unique(c(0, vapply(S, function(s) s$needsGxE, numeric(1)))))

message(sprintf("group = %s | scenarios = %d | replicates = %d | GxE founder builds/rep = %d",
                if (is.null(group)) "all" else group,
                length(S), nRep, length(gxeLevels)))

first <- TRUE
for (r in seq_len(nRep)) {
  seed <- seed0 + r
  repT0 <- Sys.time()

  ## founder objects for this replicate, one per required G x E level
  founders <- list()
  for (gg in gxeLevels) {
    pf <- base; pf$varGxERatio <- gg
    founders[[as.character(gg)]] <- makeFounders(pf, seed = seed)
  }

  for (s in S) {
    f <- founders[[as.character(s$needsGxE)]]
    p <- s$mod(base)
    p$varGxERatio <- s$needsGxE

    donors <- if (s$strategy == "CB") drawDonors(f, p, s$donorImprove) else NULL

    res <- try(runProgram(s$strategy, f, p, donors = donors, seed = seed),
               silent = TRUE)
    if (inherits(res, "try-error")) {
      message("FAILED: ", s$name, " rep ", r, " : ", as.character(res))
      next
    }

    res$scenario     <- s$name
    res$axis         <- s$axis
    res$level        <- s$level
    res$rep          <- r
    res$seed         <- seed
    res$donorImprove <- s$donorImprove
    res$group        <- s$group
    res$donorMeanObs <- if (is.null(donors)) NA_real_ else donors$realisedMean
    res$donorGapObs  <- if (is.null(donors)) NA_real_ else donors$realisedGap
    res$burnInUsed   <- f$burnInUsed
    res$eliteVarObs  <- f$eliteVarObs
    res$h2           <- p$h2
    res$cycleYears   <- p$cycleYears
    res$nBridge      <- p$nBridge
    res$bridgeBar    <- p$bridgeBar
    res$bridgeLag    <- p$bridgeLag
    res$varGxERatio  <- p$varGxERatio

    write.table(res, outfile, sep = ",", row.names = FALSE,
                col.names = first, append = !first)
    first <- FALSE
  }

  message(sprintf("rep %d/%d done in %.1f min", r, nRep,
                  as.numeric(difftime(Sys.time(), repT0, units = "mins"))))
}

message("finished -> ", outfile)
