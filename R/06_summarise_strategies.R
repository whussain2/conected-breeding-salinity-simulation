################################################################################
## 06_summarise_strategies.R -- collapse the additional-strategy runs and the
## baseline GS-RS / CB runs into the per-strategy summary used in Table S5.
##
## Usage: Rscript R/06_summarise_strategies.R
################################################################################

options(stringsAsFactors = FALSE)

grid <- do.call(rbind, lapply(Sys.glob("out/raw_[AB].csv"), function(f) {
  x <- read.csv(f); x$src <- f; x
}))
grid$repid <- paste(grid$src, grid$rep, sep = "_")
grid <- grid[grid$scenario %in% c("GSRS_base", "CB_base"), ]
grid$strategy2 <- ifelse(grid$scenario == "CB_base", "CB", "GSRS")

extra <- NULL
ef <- Sys.glob("out/raw_strategies*.csv")
if (length(ef)) {
  extra <- do.call(rbind, lapply(ef, function(f) {
    x <- read.csv(f); x$src <- f; x
  }))
  extra$repid <- paste(extra$src, extra$rep, sep = "_")
  extra$strategy2 <- extra$scenario
}

keep <- c("repid", "strategy2", "year", "gain", "varGpct", "nGeno", "nPheno")
d <- rbind(grid[, keep], if (!is.null(extra)) extra[, keep] else NULL)

## final-year row per strategy x replicate
fin <- do.call(rbind, lapply(split(d, list(d$strategy2, d$repid), drop = TRUE),
                             function(g) g[which.max(g$year), ]))

agg <- do.call(rbind, lapply(split(fin, fin$strategy2), function(g) {
  data.frame(strategy = g$strategy2[1],
             n        = nrow(g),
             gain     = sprintf("%.2f (SD %.2f)", mean(g$gain), stats::sd(g$gain)),
             gain_num = mean(g$gain),
             varpct   = sprintf("%.0f (SD %.0f)", mean(g$varGpct), stats::sd(g$varGpct)),
             varpct_num = mean(g$varGpct),
             nGeno    = round(mean(g$nGeno)),
             nPheno   = round(mean(g$nPheno)))
}))

write.csv(agg, "out/strategy_summary.csv", row.names = FALSE)
print(agg)
cat("\nwrote out/strategy_summary.csv\n")
