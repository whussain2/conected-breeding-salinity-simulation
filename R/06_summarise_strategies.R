################################################################################
## 06_summarise_strategies.R -- collapse the additional-strategy runs and the
## baseline GS-RS / CB runs into the per-strategy summary used in Table S5.
##
## Usage: Rscript R/06_summarise_strategies.R
################################################################################

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

## Guard against being pointed at the superseded first-round output, whose
## files predate the corrected bridging model and lack the genic-variance
## column.  Failing loudly here is much safer than silently producing figures
## from the wrong run.
.cb_check_output <- function(outdir) {
  f <- Sys.glob(file.path(outdir, "raw_*.csv"))
  f <- f[!grepl("strategies", f)]
  if (!length(f))
    stop("No simulation output in ", outdir, "\n",
         "  Expected files such as raw_coreA.csv. Run R/03_run.R first, or point\n",
         "  CB_SIM_ROOT at the bundle whose out/ folder contains them.",
         call. = FALSE)
  hdr <- names(utils::read.csv(f[1], nrows = 1))
  missing <- setdiff(c("scenario", "seed", "year", "gain", "varGenicPct"), hdr)
  if (length(missing))
    stop("The output in ", outdir, " is from the superseded first-round run.\n",
         "  Missing column(s): ", paste(missing, collapse = ", "), "\n",
         "  Those files predate the corrected bridging model and cannot be used\n",
         "  with this version of the scripts. Use the out/ folder that shipped\n",
         "  with this copy of the code, or re-run R/03_run.R.", call. = FALSE)
  invisible(TRUE)
}
.cb_check_output(OUTDIR)

options(stringsAsFactors = FALSE)

grid <- do.call(rbind, lapply(Sys.glob(file.path(OUTDIR, "raw_core*.csv")), function(f) {
  x <- read.csv(f); x$src <- f; x
}))
grid$repid <- paste0("seed", grid$seed)
grid <- grid[grid$scenario %in% c("GSRS_base", "CB_base"), ]
grid$strategy2 <- ifelse(grid$scenario == "CB_base", "CB", "GSRS")

extra <- NULL
ef <- Sys.glob(file.path(OUTDIR, "raw_strategies*.csv"))
if (length(ef)) {
  extra <- do.call(rbind, lapply(ef, function(f) {
    x <- read.csv(f); x$src <- f; x
  }))
  extra$repid <- paste0("seed", extra$seed)
  extra$strategy2 <- extra$scenario
}

keep <- c("repid", "strategy2", "year", "gain", "varGpct", "varGenicPct",
          "nGeno", "nPheno")
d <- rbind(grid[, keep], if (!is.null(extra)) extra[, keep] else NULL)

## final-year row per strategy x replicate
fin <- do.call(rbind, lapply(split(d, list(d$strategy2, d$repid), drop = TRUE),
                             function(g) g[which.max(g$year), ]))

agg <- do.call(rbind, lapply(split(fin, fin$strategy2), function(g) {
  data.frame(strategy = g$strategy2[1],
             n        = nrow(g),
             gain     = sprintf("%.2f (SD %.2f)", mean(g$gain), stats::sd(g$gain)),
             gain_num = mean(g$gain),
             varpct   = sprintf("%.1f (SD %.1f)", mean(g$varGenicPct),
                                stats::sd(g$varGenicPct)),
             varpct_num = mean(g$varGenicPct),
             vartot   = sprintf("%.0f (SD %.0f)", mean(g$varGpct), stats::sd(g$varGpct)),
             nGeno    = round(mean(g$nGeno)),
             nPheno   = round(mean(g$nPheno)))
}))

write.csv(agg, file.path(OUTDIR, "strategy_summary.csv"), row.names = FALSE)
print(agg)
cat("\nwrote", file.path(OUTDIR, "strategy_summary.csv"), "\n")
