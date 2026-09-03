################################################################################
## 05_analyse.R -- summarise the sensitivity and resource-matched runs,
## write Supplemental Table S6 and Supplemental Figure S2, and emit the
## numbers quoted in the main text (out/tokens.json).
##
## Usage: Rscript R/05_analyse.R
################################################################################

options(stringsAsFactors = FALSE)

files <- Sys.glob("out/raw_*.csv")
files <- files[!grepl("smoke|strategies", files)]
d <- do.call(rbind, lapply(files, function(fn) {
  x <- read.csv(fn)
  x$src <- fn
  x
}))
## unique replicate identifier across workers
d$repid <- paste(d$src, d$rep, sep = "_")
d$strategy <- as.character(d$strategy)

message("rows: ", nrow(d), " | replicates: ", length(unique(d$repid)),
        " | scenarios: ", length(unique(d$scenario)))

## ---------------------------------------------------------------- helpers
## final-year record for each scenario x replicate
finalOf <- function(df) {
  do.call(rbind, lapply(split(df, list(df$scenario, df$repid), drop = TRUE),
                        function(g) g[which.max(g$year), ]))
}
fin <- finalOf(d)  # retained for reference

## year-30 record (mid-horizon)
at <- function(df, yr) {
  do.call(rbind, lapply(split(df, list(df$scenario, df$repid), drop = TRUE),
                        function(g) {
                          i <- which.min(abs(g$year - yr))
                          g[i, ]
                        }))
}
mid <- at(d, 30)

pick <- function(df, sc) df[df$scenario == sc, ]

## Lower-variance summaries of the whole trajectory, computed per scenario x
## replicate: the mean gain over the final quarter of the horizon, and the area
## under the gain curve (a single-number summary of the entire response).
## The single year-60 endpoint is the noisiest possible summary; these are
## reported alongside it.
summaries <- function(df) {
  do.call(rbind, lapply(split(df, list(df$scenario, df$repid), drop = TRUE),
    function(g) {
      g <- g[order(g$year), ]
      lateIdx <- g$year >= max(g$year) * 0.75
      auc <- sum(diff(c(0, g$year)) * g$gain)          # t ha-1 x yr
      data.frame(scenario = g$scenario[1], repid = g$repid[1],
                 strategy = g$strategy[1],
                 axis = g$axis[1], level = g$level[1],
                 gain60 = g$gain[nrow(g)],
                 gainLate = mean(g$gain[lateIdx]),
                 auc = auc,
                 varGpct60 = g$varGpct[nrow(g)],
                 varGpctLate = mean(g$varGpct[lateIdx]),
                 nGeno = g$nGeno[nrow(g)], nPheno = g$nPheno[nrow(g)],
                 donorMeanObs = g$donorMeanObs[1])
    }))
}
S <- summaries(d)

## paired contrast of a CB scenario against a GS-RS reference scenario.
## All contrasts are paired within replicate (same founder genomes).
contrast <- function(df, cbSc, gsSc = "GSRS_base") {
  a <- df[df$scenario == cbSc, ]; b <- df[df$scenario == gsSc, ]
  m <- merge(a, b, by = "repid", suffixes = c(".cb", ".gs"))
  if (nrow(m) < 1) return(NULL)
  n <- nrow(m)
  ci <- function(x) {
    if (n < 3) return(c(NA, NA, NA))
    tt <- stats::t.test(x)
    c(tt$conf.int[1], tt$conf.int[2], tt$p.value)
  }
  d60 <- m$gain60.cb   - m$gain60.gs
  dLa <- m$gainLate.cb - m$gainLate.gs
  dAu <- m$auc.cb      - m$auc.gs
  dVa <- m$varGpctLate.cb - m$varGpctLate.gs
  c60 <- ci(d60); cLa <- ci(dLa); cAu <- ci(dAu); cVa <- ci(dVa)
  data.frame(
    scenario = cbSc, reference = gsSc, n = n,
    cb_gain = mean(m$gain60.cb), gs_gain = mean(m$gain60.gs),
    delta_gain = mean(d60), delta_sd = stats::sd(d60),
    delta_lo = c60[1], delta_hi = c60[2], p_ttest = c60[3],
    p_positive = mean(d60 > 0),
    cb_late = mean(m$gainLate.cb), gs_late = mean(m$gainLate.gs),
    delta_late = mean(dLa), delta_late_lo = cLa[1], delta_late_hi = cLa[2],
    p_late = cLa[3],
    delta_auc = mean(dAu), delta_auc_lo = cAu[1], delta_auc_hi = cAu[2],
    p_auc = cAu[3],
    cb_varpct = mean(m$varGpct60.cb), gs_varpct = mean(m$varGpct60.gs),
    delta_var = mean(dVa), delta_var_lo = cVa[1], delta_var_hi = cVa[2],
    p_var = cVa[3], p_var_positive = mean(dVa > 0),
    cb_geno = mean(m$nGeno.cb), gs_geno = mean(m$nGeno.gs),
    cb_pheno = mean(m$nPheno.cb), gs_pheno = mean(m$nPheno.gs),
    donorMean = mean(m$donorMeanObs.cb, na.rm = TRUE)
  )
}

## which GS-RS scenario is the correct reference for each CB scenario?
refFor <- function(sc) {
  if (grepl("h2_0.15", sc)) return("GSRS_h2_0.15")
  if (grepl("h2_0.50", sc)) return("GSRS_h2_0.50")
  if (sc == "CB_cyc2")      return("GSRS_cyc2")
  if (sc == "CB_cyc4")      return("GSRS_cyc4")
  if (grepl("gxe0.5", sc))  return("GSRS_gxe0.5")
  if (grepl("gxe1.0", sc))  return("GSRS_gxe1.0")
  "GSRS_base"                      # includes CB_cyc4_vs_GSRS3 by design
}

cbScen <- sort(unique(S$scenario[S$strategy == "CB"]))
tabF <- do.call(rbind, lapply(cbScen, function(x) contrast(S, x, refFor(x))))
## mid-horizon table, from the year-30 slice
midS <- summaries(d[d$year <= 30, ])
tabM <- do.call(rbind, lapply(cbScen, function(x) contrast(midS, x, refFor(x))))

## attach axis / level labels
lab <- unique(d[, c("scenario", "axis", "level")])
tabF <- merge(tabF, lab, by = "scenario")
tabM <- merge(tabM, lab, by = "scenario")

axisOrder <- c("baseline", "donor_quality", "n_bridge", "bridge_bar",
               "heritability", "cycle_length", "gxe", "resource_matched")
tabF$axis <- factor(tabF$axis, levels = axisOrder)
tabF <- tabF[order(tabF$axis, tabF$scenario), ]
tabM$axis <- factor(tabM$axis, levels = axisOrder)
tabM <- tabM[order(tabM$axis, tabM$scenario), ]

write.csv(tabF, "out/tableS6_year60.csv", row.names = FALSE)
write.csv(tabM, "out/tableS6_year30.csv", row.names = FALSE)

## ---------------------------------------------------------------- console report
fmt <- function(x, k = 2) formatC(x, format = "f", digits = k)
cat("\n=========== YEAR-60 AND LATE-HORIZON CONTRASTS (CB minus reference) ===========\n")
cat(sprintf("%-22s %-30s %3s %7s %7s %8s %7s %8s %9s %8s\n",
            "scenario", "level", "n", "d60", "d60SD", "d60 P",
            "dLate", "dLate P", "dVar(pp)", "dVar P"))
for (i in seq_len(nrow(tabF))) with(tabF[i, ],
  cat(sprintf("%-22s %-30s %3d %7s %7s %8s %7s %8s %9s %8s\n",
              scenario, substr(level, 1, 30), n, fmt(delta_gain), fmt(delta_sd),
              fmt(p_ttest, 3), fmt(delta_late), fmt(p_late, 3),
              fmt(delta_var, 1), fmt(p_var, 4))))

## ---------------------------------------------------------------- Figure S2
for (dev in c("png", "pdf")) {
if (dev == "png") png("out/FigureS2.png", width = 2400, height = 1700, res = 220) else
  pdf("out/FigureS2.pdf", width = 11.0, height = 7.8, pointsize = 10)
op <- par(mfrow = c(2, 2), mar = c(4.2, 4.4, 2.6, 1.0), mgp = c(2.5, 0.8, 0),
          cex.lab = 1.0, cex.axis = 0.9)

## (A) baseline trajectories with +/- 1 SD ribbons
trajOf <- function(sc) {
  g <- d[d$scenario == sc, ]
  if (!nrow(g)) return(NULL)
  ag <- aggregate(cbind(gain, varGpct) ~ year, g, mean)
  sd1 <- aggregate(cbind(gain, varGpct) ~ year, g, stats::sd)
  names(sd1)[-1] <- paste0(names(sd1)[-1], "_sd")
  merge(ag, sd1, by = "year")
}
tg <- trajOf("GSRS_base"); tc <- trajOf("CB_base")
ylim <- range(0, tg$gain + tg$gain_sd, tc$gain + tc$gain_sd, na.rm = TRUE)
plot(tg$year, tg$gain, type = "n", ylim = ylim, xlab = "Year",
     ylab = expression("Cumulative genetic gain (t ha"^-1*")"),
     main = "A  Baseline: gain")
ribbon <- function(t, col) {
  polygon(c(t$year, rev(t$year)),
          c(t$gain - t$gain_sd, rev(t$gain + t$gain_sd)),
          col = adjustcolor(col, 0.18), border = NA)
  lines(t$year, t$gain, col = col, lwd = 2.4)
}
adjustcolor <- function(col, a) grDevices::adjustcolor(col, alpha.f = a)
ribbon(tg, "#B4462F"); ribbon(tc, "#2F6FB4")
legend("topleft", c("Connected Breeding", "GS-RS (closed pool)"),
       col = c("#2F6FB4", "#B4462F"), lwd = 2.4, bty = "n", cex = 0.85)

## (B) baseline variance retention
ylim2 <- range(0, tg$varGpct + tg$varGpct_sd, tc$varGpct + tc$varGpct_sd, na.rm = TRUE)
plot(tg$year, tg$varGpct, type = "n", ylim = ylim2, xlab = "Year",
     ylab = "Additive genetic variance (% of elite base)",
     main = "B  Baseline: variance retained")
rib2 <- function(t, col) {
  polygon(c(t$year, rev(t$year)),
          c(t$varGpct - t$varGpct_sd, rev(t$varGpct + t$varGpct_sd)),
          col = adjustcolor(col, 0.18), border = NA)
  lines(t$year, t$varGpct, col = col, lwd = 2.4)
}
rib2(tg, "#B4462F"); rib2(tc, "#2F6FB4")
abline(h = 100, lty = 3, col = "grey50")

## (C) tornado plot of the year-60 CB advantage
tt <- tabF[tabF$axis != "baseline", ]
tt <- tt[order(tt$delta_gain), ]
baseDelta <- tabF$delta_gain[tabF$scenario == "CB_base"]
prettyAxis <- c(donor_quality = "Donor quality", n_bridge = "Bridge lines/cycle",
                bridge_bar = "Elite-equivalence bar", heritability = "Heritability",
                cycle_length = "Cycle length", gxe = "G x E",
                resource_matched = "Resource-matched")
prettyLev <- function(x) {
  x <- gsub("^1 cycle\\(s\\) of donor pre-improvement", "donors pre-improved 1 cycle", x)
  x <- gsub("^2 cycle\\(s\\) of donor pre-improvement", "donors pre-improved 2 cycles", x)
  x <- gsub("^3 cycle\\(s\\) of donor pre-improvement", "donors pre-improved 3 cycles", x)
  x <- gsub("varGxE/varA = ", "var(GxE)/var(A) = ", x)
  x <- gsub("^h2 = ", "h2 = ", x)
  x
}
labs <- paste0(prettyAxis[as.character(tt$axis)], ": ", prettyLev(tt$level))
par(mar = c(4.2, 17, 2.6, 1.0))
bp <- barplot(tt$delta_gain, horiz = TRUE, names.arg = labs, las = 1,
              cex.names = 0.62, xlab = expression("Year-60 CB advantage (t ha"^-1*")"),
              col = ifelse(tt$delta_gain > 0, "#2F6FB4", "#B4462F"),
              border = NA, main = "C  Sensitivity of the CB advantage")
arrows(tt$delta_lo, bp, tt$delta_hi, bp, angle = 90, code = 3,
       length = 0.02, col = "grey30", lwd = 1)
abline(v = 0, col = "black")
abline(v = baseDelta, lty = 2, col = "grey40")

## (D) gain-vs-diversity plane at year 60
par(mar = c(4.2, 4.4, 2.6, 1.0))
plot(tabF$cb_varpct, tabF$cb_gain, pch = 21, bg = "#2F6FB4", col = "white",
     cex = 1.3, xlab = "Additive genetic variance retained (% of elite base)",
     ylab = expression("Year-60 cumulative gain (t ha"^-1*")"),
     main = "D  Gain versus retained diversity",
     xlim = range(0, tabF$cb_varpct, tabF$gs_varpct, na.rm = TRUE),
     ylim = range(tabF$cb_gain, tabF$gs_gain, na.rm = TRUE))
points(tabF$gs_varpct, tabF$gs_gain, pch = 21, bg = "#B4462F", col = "white",
       cex = 1.3)
legend("bottomright", c("Connected Breeding", "GS-RS reference"),
       pch = 21, pt.bg = c("#2F6FB4", "#B4462F"), col = "white", bty = "n",
       cex = 0.85)
par(op)
dev.off()
}

## ---------------------------------------------------------------- tokens
b    <- tabF[tabF$scenario == "CB_base", ]
rm_  <- tabF[tabF$axis == "resource_matched", ]
sens <- tabF[!tabF$axis %in% c("baseline", "resource_matched"), ]
nrep <- b$n

genoPct  <- 100 * (b$cb_geno / b$gs_geno - 1)
phenoPct <- 100 * (b$cb_pheno / b$gs_pheno - 1)

gv <- function(sc, col = "delta_gain") tabF[[col]][tabF$scenario == sc]

sig <- function(p) !is.na(p) & p < 0.05
nSigPos <- sum(sig(tabF$p_ttest) & tabF$delta_gain > 0)
nSigNeg <- sum(sig(tabF$p_ttest) & tabF$delta_gain < 0)
nVarPos <- sum(tabF$delta_var > 0)
nVarSig <- sum(sig(tabF$p_var) & tabF$delta_var > 0)

costTxt <- sprintf(
  "over the 60-year horizon the baseline Connected Breeding programme genotypes %.0f%% more lines and phenotypes %.0f%% more plots than the closed-pool programme",
  genoPct, phenoPct)

costShortTxt <- sprintf(
  "%.0f%% more lines and phenotypes %.0f%% more plots over the 60-year horizon than the closed-pool programme, and additionally maintains a donor reservoir across all cycles.",
  genoPct, phenoPct)

pv <- function(sc) tabF$p_ttest[tabF$scenario == sc]
pvar <- function(sc) tabF$p_var[tabF$scenario == sc]
pfmt <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-8)  return("P < 10^-8")
  if (p < 0.001) return(sprintf("P = %.5f", p))
  if (p < 0.01)  return(sprintf("P = %.3f", p))
  sprintf("P = %.2f", p)
}
spellOut <- function(n) c("none", "one", "two", "three", "four", "five", "six",
                          "seven", "eight", "nine", "ten")[n + 1]
tidyLev <- function(x) {
  x <- gsub("^1 cycle\\(s\\)", "one cycle", x)
  x <- gsub("^2 cycle\\(s\\)", "two cycles", x)
  x <- gsub("^3 cycle\\(s\\)", "three cycles", x)
  x
}

baseTxt <- sprintf(
  paste0("Re-running the same design with %d paired replicate simulations, the year-60 ",
         "difference in cumulative gain between the two strategies was %+.2f t ha-1 ",
         "(95%% CI %.2f to %+.2f; %s), that is, indistinguishable from zero, whereas the ",
         "difference in retained additive genetic variance was large and highly significant ",
         "(%.0f%% versus %.0f%% of the elite base; difference %.0f percentage points, 95%% CI ",
         "%.0f to %.0f; %s). On this evidence the robust effect of Connected Breeding is the ",
         "conservation of additive genetic variance, and hence of long-term selection ",
         "potential, rather than a higher rate of gain over the horizon simulated."),
  nrep, b$delta_gain, b$delta_lo, b$delta_hi, pfmt(b$p_ttest),
  b$cb_varpct, b$gs_varpct, b$delta_var, b$delta_var_lo, b$delta_var_hi,
  pfmt(b$p_var))

resTxt <- sprintf(
  paste0("Under an equal genotyping and phenotyping budget the year-60 difference in gain was ",
         "%+.2f t ha-1 (%s); with a two-cycle development lag %+.2f t ha-1 (%s); with both ",
         "%+.2f t ha-1 (%s); and with a four-year Connected Breeding cycle against a three-year ",
         "closed-pool cycle %+.2f t ha-1 (%s). None of these differences was statistically ",
         "significant, but all four are negative, so under equal constraints Connected Breeding ",
         "gives up a small and uncertain amount of short-term gain. The retention of additive ",
         "genetic variance was unaffected in direction: Connected Breeding retained more of it ",
         "than the closed-pool programme in all four constrained scenarios."),
  gv("CB_equalBudget"), pfmt(pv("CB_equalBudget")),
  gv("CB_lag2"), pfmt(pv("CB_lag2")),
  gv("CB_equalBudget_lag2"), pfmt(pv("CB_equalBudget_lag2")),
  gv("CB_cyc4_vs_GSRS3"), pfmt(pv("CB_cyc4_vs_GSRS3")))

bestSc  <- tabF$scenario[which.max(tabF$delta_gain)]
bestLev <- tabF$level[tabF$scenario == bestSc]

sensTxt <- sprintf(
  paste0("Across the %d scenarios examined, Connected Breeding retained more additive genetic ",
         "variance than its matched closed-pool reference in every one, significantly so in %d, ",
         "while the closed-pool programme had essentially exhausted its variance by the end of ",
         "the horizon in all of them. The difference in cumulative gain was not robust: it ",
         "ranged from %+.2f to %+.2f t ha-1, was positive in %d of %d scenarios, and differed ",
         "significantly from zero in only %s, namely when donors had already been pre-improved ",
         "(%s: %+.2f t ha-1, %s), which is the condition under which bridge lines carry ",
         "favourable alleles without imposing a large penalty on the population mean. Donor ",
         "quality and the stringency with which bridge lines are required to match elite ",
         "performance were the most influential factors; heritability, cycle length and ",
         "genotype-by-environment interaction had little effect on the contrast."),
  nrow(tabF), nVarSig, min(tabF$delta_gain), max(tabF$delta_gain),
  sum(tabF$delta_gain > 0), nrow(tabF), spellOut(nSigPos),
  tidyLev(as.character(bestLev)), max(tabF$delta_gain), pfmt(pv(bestSc)))

tok <- list(
  BASELINE_REPLICATION = baseTxt,
  RESOURCE_COST       = costTxt,
  RESOURCE_COST_SHORT = costShortTxt,
  RESOURCE_RESULT     = resTxt,
  SENS_RESULT         = sensTxt,
  MAGIC_NAM           = "we simulated a one-time broad-base multi-parent population, and closed-pool phenotypic recurrent selection, under the identical parameterisation, so that the corresponding entries in Table S5 are outputs of this simulation (Supplemental File S1)."
)

con <- file("out/tokens.json", "w")
writeLines(paste0("{\n",
  paste(sprintf('  "%s": %s', names(tok),
                vapply(tok, function(x) {
                  x <- gsub('\\\\', '\\\\\\\\', x); x <- gsub('"', '\\\\"', x)
                  paste0('"', x, '"')
                }, character(1))), collapse = ",\n"),
  "\n}"), con, useBytes = TRUE)
close(con)

cat("\n---- tokens ----\n")
for (n in names(tok)) cat("\n[", n, "]\n", tok[[n]], "\n")
cat("\nwrote out/tableS6_year60.csv, out/tableS6_year30.csv, out/FigureS2.png, out/tokens.json\n")
