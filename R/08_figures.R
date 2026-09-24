################################################################################
## 08_figures.R -- every simulation-derived figure panel in the paper.
##
## Writes, in both PDF (vector, for submission) and PNG (300 dpi preview):
##   figures/Figure_4C_gain.*              cumulative genetic gain, CB vs GS-RS
##   figures/Figure_4D_diversity.*         genic variance retained, CB vs GS-RS
##   figures/Figure_S2_A_inflow_gain.*     gain across the donor-inflow ladder
##   figures/Figure_S2_B_inflow_diversity.*  diversity across the ladder
##   figures/Figure_S2_C_late_slope.*      late-horizon rate of gain
##   figures/Figure_S2_D_sensitivity.*     year-60 contrasts, all scenarios
##   figures/Figure_4CD_combined.*         C and D side by side
##   figures/Figure_S2_combined.*          the four S2 panels on one sheet
##
## and the exact plotted values as CSV alongside each figure, so every panel is
## reproducible from data without re-running the simulation.
##
## Usage: Rscript R/08_figures.R
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
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------- palette
## Validated categorical/sequential palette (see README).  GS-RS takes the
## contrasting categorical hue; the Connected Breeding inflow ladder is an
## ORDERED variable, so it takes steps of a single sequential blue ramp.
COL <- list(
  gsrs   = "#eb6834",
  cb1    = "#86b6ef",
  cb2    = "#3987e5",
  cb4    = "#256abf",
  cb8    = "#104281",
  ink    = "#0b0b0b",
  ink2   = "#52514e",
  muted  = "#8a8983",
  grid   = "#e6e5e0",
  neg    = "#c0442a",
  pos    = "#2a78d6"
)
alpha <- function(col, a) grDevices::adjustcolor(col, alpha.f = a)

## ---------------------------------------------------------------- data
files <- Sys.glob(file.path(OUTDIR, "raw_*.csv"))
files <- files[!grepl("strategies", files)]
d <- do.call(rbind, lapply(files, function(fn) {
  x <- read.csv(fn)
  keep <- c("scenario", "axis", "level", "strategy", "rep", "seed", "year",
            "gain", "varG", "varGpct", "varGenic", "varGenicPct",
            "nGeno", "nPheno", "nElig", "nBridgeIn")
  x <- x[, intersect(keep, names(x))]
  x$src <- basename(fn)
  x
}))
## keyed by seed so scenarios from different batches pair correctly
d$repid <- paste0("seed", d$seed)
message("rows ", nrow(d), " | replicates ", length(unique(d$repid)),
        " | scenarios ", length(unique(d$scenario)))

traj <- function(sc, value) {
  g <- d[d$scenario == sc, ]
  if (!nrow(g)) return(NULL)
  m  <- tapply(g[[value]], g$year, mean)
  s  <- tapply(g[[value]], g$year, stats::sd)
  n  <- tapply(g[[value]], g$year, length)
  data.frame(year = as.numeric(names(m)), mean = as.numeric(m),
             sd = as.numeric(s), n = as.numeric(n))
}

## per-replicate late-horizon slope of cumulative gain
slopes <- function(sc, from = 45) {
  g <- d[d$scenario == sc, ]
  if (!nrow(g)) return(numeric(0))
  vapply(split(g, g$repid), function(x) {
    x <- x[x$year >= from, ]
    if (nrow(x) < 3) return(NA_real_)
    unname(stats::coef(stats::lm(gain ~ year, x))[2])
  }, numeric(1))
}

final <- function(sc, value = "gain") {
  g <- d[d$scenario == sc, ]
  if (!nrow(g)) return(numeric(0))
  vapply(split(g, g$repid), function(x) x[[value]][which.max(x$year)], numeric(1))
}

## ---------------------------------------------------------------- devices
open_dev <- function(stem, w, h) {
  list(
    function() pdf(file.path(FIGDIR, paste0(stem, ".pdf")), width = w,
                   height = h, pointsize = 9),
    function() png(file.path(FIGDIR, paste0(stem, ".png")),
                   width = w * 300, height = h * 300, res = 300)
  )
}
render <- function(stem, w, h, fn) {
  for (mk in open_dev(stem, w, h)) { mk(); fn(); dev.off() }
  message("wrote ", file.path(basename(FIGDIR), paste0(stem, ".pdf/.png")))
}

axis_style <- function() {
  box(col = COL$muted, lwd = 0.8)
}
gridlines <- function(h = TRUE) {
  if (h) abline(h = axTicks(2), col = COL$grid, lwd = 0.7)
}

ribbon <- function(t, col, lwd = 2) {
  ok <- !is.na(t$sd)
  lg <- par("ylog")
  floorv <- if (lg) 10^par("usr")[3] else -Inf
  lower <- pmax((t$mean - t$sd)[ok], floorv)
  polygon(c(t$year[ok], rev(t$year[ok])),
          c(lower, rev((t$mean + t$sd)[ok])),
          col = alpha(col, 0.16), border = NA)
  lines(t$year, t$mean, col = col, lwd = lwd)
}

## place a set of end-of-line labels, nudging them apart so they never collide
## Nudge end-of-line labels apart so they never collide.  Separation is applied
## in plotting coordinates, which on a log axis means log10 units, and the
## minimum gap defaults to one line height.
endlabels <- function(items, minsep = NULL, cex = 0.72) {
  lg <- par("ylog")
  usr <- par("usr")[3:4]
  if (is.null(minsep)) minsep <- diff(usr) * 0.055 * cex / 0.72
  raw <- vapply(items, function(x) x$t$mean[nrow(x$t)], 0)
  ys  <- if (lg) log10(raw) else raw
  o   <- order(ys)
  ys2 <- ys
  for (k in seq_along(o)[-1]) {
    i <- o[k]; j <- o[k - 1]
    if (ys2[i] - ys2[j] < minsep) ys2[i] <- ys2[j] + minsep
  }
  for (i in seq_along(items))
    text(items[[i]]$t$year[nrow(items[[i]]$t)],
         if (lg) 10^ys2[i] else ys2[i], items[[i]]$lab,
         col = items[[i]]$col, pos = 4, offset = 0.3, cex = cex, xpd = NA,
         font = 2)
}

## ================================================================ Figure 4C
gs <- traj("GSRS_base", "gain"); cb <- traj("CB_base", "gain")
stopifnot(!is.null(gs), !is.null(cb))

panel_4C <- function(withLadder = FALSE) {
  ys <- c(gs$mean + gs$sd, cb$mean + cb$sd, 0)
  plot(NA, xlim = c(0, 68), ylim = c(0, max(ys, na.rm = TRUE) * 1.05),
       xlab = "Year", ylab = expression("Cumulative genetic gain (t ha"^-1*")"),
       xaxs = "i", yaxs = "i", axes = FALSE)
  axis(1, col = COL$muted, col.axis = COL$ink2, at = seq(0, 60, 10))
  axis(2, col = COL$muted, col.axis = COL$ink2, las = 1)
  gridlines(); axis_style()
  ribbon(gs, COL$gsrs); ribbon(cb, COL$cb2)
  endlabels(list(list(t = cb, lab = "Connected\nBreeding", col = COL$cb2),
                 list(t = gs, lab = "GS-RS", col = COL$gsrs)),
            )
  legend("topleft", c("Connected Breeding", "GS-RS (closed pool)"),
         col = c(COL$cb2, COL$gsrs), lwd = 2, bty = "n", cex = 0.78,
         text.col = COL$ink2, seg.len = 1.4)
}
render("Figure_4C_gain", 4.6, 3.5, function() {
  par(mar = c(3.6, 4.0, 1.0, 5.6), mgp = c(2.3, 0.6, 0), tcl = -0.25,
      col.lab = COL$ink, family = "sans")
  panel_4C()
})
write.csv(rbind(cbind(strategy = "GS-RS", gs), cbind(strategy = "CB", cb)),
          file.path(FIGDIR, "Figure_4C_gain_data.csv"), row.names = FALSE)

## ================================================================ Figure 4D
gsv <- traj("GSRS_base", "varGenicPct"); cbv <- traj("CB_base", "varGenicPct")

## A logarithmic axis is used here deliberately: the quantity falls through more
## than one order of magnitude, and on a linear axis the difference that matters
## -- the several-fold gap in the diversity still available at the end of the
## horizon -- is compressed into invisibility near zero.
panel_4D <- function() {
  ticks <- c(1, 2, 5, 10, 20, 50, 100)
  lo <- min(c(gsv$mean, cbv$mean), na.rm = TRUE)
  plot(NA, xlim = c(0, 68), ylim = c(max(1, lo * 0.7), 115), log = "y",
       xlab = "Year", ylab = "Genic variance retained (% of elite base)",
       xaxs = "i", axes = FALSE)
  axis(1, col = COL$muted, col.axis = COL$ink2, at = seq(0, 60, 10))
  axis(2, at = ticks, labels = ticks, col = COL$muted, col.axis = COL$ink2,
       las = 1)
  abline(h = ticks, col = COL$grid, lwd = 0.7); axis_style()
  ribbon(gsv, COL$gsrs); ribbon(cbv, COL$cb2)
  endlabels(list(list(t = cbv, lab = "Connected\nBreeding", col = COL$cb2),
                 list(t = gsv, lab = "GS-RS", col = COL$gsrs)))
  legend("bottomleft", c("Connected Breeding", "GS-RS (closed pool)"),
         col = c(COL$cb2, COL$gsrs), lwd = 2, bty = "n", cex = 0.78,
         text.col = COL$ink2, seg.len = 1.4)
}
render("Figure_4D_diversity", 4.6, 3.5, function() {
  par(mar = c(3.6, 4.2, 1.0, 5.6), mgp = c(2.5, 0.6, 0), tcl = -0.25,
      col.lab = COL$ink, family = "sans")
  panel_4D()
})
write.csv(rbind(cbind(strategy = "GS-RS", gsv), cbind(strategy = "CB", cbv)),
          file.path(FIGDIR, "Figure_4D_diversity_data.csv"), row.names = FALSE)

render("Figure_4CD_combined", 9.2, 3.6, function() {
  par(mfrow = c(1, 2), mar = c(3.6, 4.2, 1.8, 5.4), mgp = c(2.4, 0.6, 0),
      tcl = -0.25, col.lab = COL$ink, family = "sans")
  panel_4C(); mtext("C", side = 3, adj = -0.16, line = 0.4, font = 2, cex = 1.1)
  panel_4D(); mtext("D", side = 3, adj = -0.16, line = 0.4, font = 2, cex = 1.1)
})

## ================================================ Figure S2 A/B: inflow ladder
LADDER <- list(
  list(sc = "GSRS_base", lab = "GS-RS (no inflow)",       col = COL$gsrs),
  list(sc = "CB_nB1",    lab = "1 bridge line / cycle",   col = COL$cb1),
  list(sc = "CB_base",   lab = "2 (as specified)",        col = COL$cb2),
  list(sc = "CB_nB4",    lab = "4",                       col = COL$cb4),
  list(sc = "CB_nB8",    lab = "8",                       col = COL$cb8)
)
LADDER <- LADDER[vapply(LADDER, function(x) !is.null(traj(x$sc, "gain")),
                        logical(1))]

ladder_panel <- function(value, ylab, legpos = "topleft", logy = FALSE) {
  ts <- lapply(LADDER, function(x) traj(x$sc, value))
  ymax <- max(vapply(ts, function(t) max(t$mean, na.rm = TRUE), 0)) * 1.08
  if (logy) {
    ticks <- c(1, 2, 5, 10, 20, 50, 100)
    lo <- min(vapply(ts, function(t) min(t$mean, na.rm = TRUE), Inf))
    plot(NA, xlim = c(0, 70), ylim = c(max(1, lo * 0.7), 115), log = "y",
         xlab = "Year", ylab = ylab, xaxs = "i", axes = FALSE)
    axis(1, col = COL$muted, col.axis = COL$ink2, at = seq(0, 60, 10))
    axis(2, at = ticks, labels = ticks, col = COL$muted, col.axis = COL$ink2,
         las = 1)
    abline(h = ticks, col = COL$grid, lwd = 0.7)
  } else {
  plot(NA, xlim = c(0, 70), ylim = c(0, ymax), xlab = "Year", ylab = ylab,
       xaxs = "i", yaxs = "i", axes = FALSE)
  axis(1, col = COL$muted, col.axis = COL$ink2, at = seq(0, 60, 10))
  axis(2, col = COL$muted, col.axis = COL$ink2, las = 1)
  gridlines()
  }
  axis_style()
  for (i in seq_along(LADDER))
    lines(ts[[i]]$year, ts[[i]]$mean, col = LADDER[[i]]$col,
          lwd = if (LADDER[[i]]$sc %in% c("GSRS_base", "CB_base")) 2.4 else 1.6)
  endlabels(lapply(seq_along(LADDER), function(i)
              list(t = ts[[i]], lab = LADDER[[i]]$lab, col = LADDER[[i]]$col)),
            cex = 0.68)
  legend(legpos, legend = vapply(LADDER, function(x) x$lab, ""),
         col = vapply(LADDER, function(x) x$col, ""), lwd = 2, bty = "n",
         cex = 0.62, text.col = COL$ink2, seg.len = 1.2,
         title = "Bridge lines admitted per cycle",
         title.col = COL$ink, title.adj = 0, title.cex = 0.62)
  invisible(ts)
}

render("Figure_S2_A_inflow_gain", 5.0, 3.6, function() {
  par(mar = c(3.6, 4.0, 1.0, 6.6), mgp = c(2.3, 0.6, 0), tcl = -0.25,
      col.lab = COL$ink, family = "sans")
  ladder_panel("gain", expression("Cumulative genetic gain (t ha"^-1*")"))
})
render("Figure_S2_B_inflow_diversity", 5.0, 3.6, function() {
  par(mar = c(3.6, 4.2, 1.0, 6.6), mgp = c(2.5, 0.6, 0), tcl = -0.25,
      col.lab = COL$ink, family = "sans")
  ladder_panel("varGenicPct", "Genic variance retained (% of elite base)",
               "bottomleft", logy = TRUE)
})
lad <- do.call(rbind, lapply(LADDER, function(x) {
  a <- traj(x$sc, "gain"); b <- traj(x$sc, "varGenicPct")
  data.frame(scenario = x$sc, label = x$lab, year = a$year,
             gain = a$mean, gain_sd = a$sd,
             genicPct = b$mean, genicPct_sd = b$sd, n = a$n)
}))
write.csv(lad, file.path(FIGDIR, "Figure_S2_AB_inflow_data.csv"), row.names = FALSE)

## ============================================ Figure S2 C: late-horizon slope
SLOPE_SET <- LADDER
slp <- lapply(SLOPE_SET, function(x) slopes(x$sc, 45))
slp_tab <- do.call(rbind, lapply(seq_along(SLOPE_SET), function(i) {
  v <- slp[[i]][!is.na(slp[[i]])]
  if (!length(v)) return(NULL)
  ci <- if (length(v) > 2) stats::t.test(v)$conf.int else c(NA, NA)
  data.frame(scenario = SLOPE_SET[[i]]$sc, label = SLOPE_SET[[i]]$lab,
             n = length(v), slope = mean(v), sd = stats::sd(v),
             lo = ci[1], hi = ci[2], col = SLOPE_SET[[i]]$col)
}))

panel_S2C <- function() {
  s <- slp_tab
  xm <- max(c(s$hi, s$slope), na.rm = TRUE) * 1.30
  bp <- barplot(s$slope, horiz = TRUE, names.arg = rep("", nrow(s)),
                col = s$col, border = NA,
                xlim = c(0, xm),
                xlab = expression("Rate of gain, years 45-60 (t ha"^-1*" yr"^-1*")"),
                axes = FALSE)
  axis(1, col = COL$muted, col.axis = COL$ink2)
  abline(v = 0, col = COL$ink, lwd = 0.9)
  suppressWarnings(arrows(s$lo, bp, s$hi, bp, angle = 90, code = 3,
                          length = 0.03, col = COL$ink2, lwd = 0.9))
  ## category labels sit in the left margin, values beyond the whisker
  text(par("usr")[1] - diff(par("usr")[1:2]) * 0.02, bp, s$label,
       adj = 1, cex = 0.7, col = COL$ink, xpd = NA)
  text(pmax(s$hi, s$slope, na.rm = TRUE), bp, sprintf("%+.4f", s$slope),
       pos = 4, offset = 0.3, cex = 0.66, col = COL$ink2, xpd = NA)
}
render("Figure_S2_C_late_slope", 6.2, 3.0, function() {
  par(mar = c(3.8, 11.5, 1.0, 4.0), mgp = c(2.4, 0.6, 0), tcl = -0.25,
      col.lab = COL$ink, family = "sans")
  panel_S2C()
})
write.csv(slp_tab[, setdiff(names(slp_tab), "col")],
          file.path(FIGDIR, "Figure_S2_C_late_slope_data.csv"), row.names = FALSE)

## ======================================= Figure S2 D: year-60 contrast, all
refFor <- function(sc) {
  if (grepl("h2_0.15", sc)) return("GSRS_h2_0.15")
  if (grepl("h2_0.50", sc)) return("GSRS_h2_0.50")
  if (sc == "CB_cyc2")      return("GSRS_cyc2")
  if (sc == "CB_cyc4")      return("GSRS_cyc4")
  if (grepl("gxe0.5", sc))  return("GSRS_gxe0.5")
  if (grepl("gxe1.0", sc))  return("GSRS_gxe1.0")
  "GSRS_base"
}
AXLAB <- c(baseline = "Baseline", inflow_rate = "Donor inflow rate",
           mechanism = "Mechanism control", bridge_bar = "Elite-equivalence bar",
           donor_quality = "Donor quality", heritability = "Heritability",
           cycle_length = "Cycle length", gxe = "G x E",
           resource_matched = "Resource-matched")

cbScen <- unique(d$scenario[d$strategy == "CB"])
contr <- do.call(rbind, lapply(cbScen, function(sc) {
  a <- final(sc); b <- final(refFor(sc))
  common <- intersect(names(a), names(b))
  if (length(common) < 2) return(NULL)
  v <- a[common] - b[common]
  ci <- if (length(v) > 2) stats::t.test(v)$conf.int else c(NA, NA)
  pv <- if (length(v) > 2) stats::t.test(v)$p.value else NA
  info <- d[d$scenario == sc, ][1, ]
  data.frame(scenario = sc, axis = info$axis, level = info$level,
             n = length(v), delta = mean(v), sd = stats::sd(v),
             lo = ci[1], hi = ci[2], p = pv)
}))
contr$lab <- paste0(AXLAB[contr$axis], ": ", contr$level)
contr <- contr[order(contr$delta), ]

panel_S2D <- function() {
  s <- contr
  xr <- range(c(s$lo, s$hi, s$delta, 0), na.rm = TRUE)
  bp <- barplot(s$delta, horiz = TRUE, names.arg = rep("", nrow(s)),
                col = ifelse(s$delta > 0, COL$pos, COL$neg), border = NA,
                xlim = xr + c(-0.05, 0.05) * diff(xr),
                xlab = expression("Year-60 difference vs closed pool (t ha"^-1*")"),
                axes = FALSE)
  axis(1, col = COL$muted, col.axis = COL$ink2)
  suppressWarnings(arrows(s$lo, bp, s$hi, bp, angle = 90, code = 3,
                          length = 0.02, col = COL$ink2, lwd = 0.8))
  abline(v = 0, col = COL$ink, lwd = 0.9)
  text(par("usr")[1] - diff(par("usr")[1:2]) * 0.015, bp, s$lab,
       adj = 1, cex = 0.58, col = COL$ink, xpd = NA)
}
render("Figure_S2_D_sensitivity", 8.2, 0.24 * nrow(contr) + 1.4, function() {
  par(mar = c(3.6, 19.0, 0.8, 1.2), mgp = c(2.3, 0.6, 0), tcl = -0.25,
      col.lab = COL$ink, family = "sans")
  panel_S2D()
})
write.csv(contr, file.path(FIGDIR, "Figure_S2_D_sensitivity_data.csv"), row.names = FALSE)

## ---------------------------------------------------------------- combined S2
render("Figure_S2_combined", 12.0, 9.0, function() {
  layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
  par(mar = c(3.8, 4.2, 2.2, 6.4), mgp = c(2.4, 0.6, 0), tcl = -0.25,
      col.lab = COL$ink, family = "sans")
  ladder_panel("gain", expression("Cumulative genetic gain (t ha"^-1*")"))
  mtext("A", side = 3, adj = -0.14, line = 0.6, font = 2, cex = 1.1)
  ladder_panel("varGenicPct", "Genic variance retained (% of elite base)",
               "bottomleft", logy = TRUE)
  mtext("B", side = 3, adj = -0.14, line = 0.6, font = 2, cex = 1.1)
  par(mar = c(3.8, 10.5, 2.2, 4.0))
  panel_S2C()
  mtext("C", side = 3, adj = -0.42, line = 0.6, font = 2, cex = 1.1)
  par(mar = c(3.8, 16.5, 2.2, 1.2))
  panel_S2D()
  mtext("D", side = 3, adj = -0.72, line = 0.6, font = 2, cex = 1.1)
})

message("\nAll figure panels and their plotted data are in ", FIGDIR)
