################################################################################
## 08_figures2.R -- Revised script for creating figures for manuscript revision, August 2024
##
## Each block below reads its own CSV and builds its own plot. Nothing is shared
## between blocks, so you can run any one of them on its own, copy it into a new
## script, or delete the ones you do not want. Colours, labels and sizes are
## written inline where they are used.
##
## Run this once at the start of the session, then run whichever block you want:
################################################################################

## needs ggplot2 once:  install.packages("ggplot2")
library(ggplot2)

## Point R at the bundle folder (the one containing figures/ and out/).
## This is the only line you need to edit.
setwd("~/Documents/IRRI_Files/12_Manuscript/Salinity_TrendsInPlantScince/Plant Communications/Revision_2nd_August19/Secon_simulation/connected-breeding-salinity-simulation/connected-breeding-salinity-simulation")

dir.create("figures2", showWarnings = FALSE)   # where the new plots are written


################################################################################
## FIGURE 4C -- cumulative genetic gain, Connected Breeding vs GS-RS
## columns: strategy, year, mean, sd, n
################################################################################

d <- read.csv("figures/Figure_4C_gain_data.csv")
head(d)

## legend row order: Connected Breeding first, GS-RS second
d$strategy <- factor(d$strategy, levels = c("CB", "GS-RS"))

png("Figure4c_MYMYYY.png", units="in", width=10, height=6, res=600)
p_4C <- ggplot(d, aes(x = year, y = mean, colour = strategy, fill = strategy)) +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +                     # dots take the line's colour
  scale_colour_manual(values = c("CB" = "#2d6d66", "GS-RS" = "#a0522d"),
                      labels = c("Connected Breeding (GS + bridge materials)", "Recurrent genomic selection (closed pool)")) +
  scale_fill_manual(values = c("CB" = "#2d6d66", "GS-RS" = "#a0522d"),
                    labels = c("Connected Breeding (GS + bridge materials)", "Recurrent genomic selection (closed pool)")) +
  labs(x = "", y = "",
       colour = NULL, fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    ## legend inside the panel, upper left.
    ## legend.position   = c(x, y): fractions of the PANEL, 0 = left/bottom,
    ##                     1 = right/top. Nudge these two numbers to move it.
    ## legend.justification = which corner of the legend box sits on that point;
    ##                     c(0, 1) = the legend's top-left corner.
    legend.position      = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.key           = element_blank(),
    legend.spacing.y     = unit(0, "pt"),
    legend.margin        = margin(2, 4, 2, 2)
  )

p_4C
dev.off()
#ggsave("figures2/Fig_4C.pdf", p_4C, width = 6, height = 4.5)


################################################################################
## FIGURE 4D -- genic variance retained, same two strategies
## columns: strategy, year, mean, sd, n
## Log y-axis: the quantity falls through more than one order of magnitude, and
## the gap that matters is the one left at the end. Delete scale_y_log10() for
## a linear axis.
################################################################################

d <- read.csv("figures/Figure_4D_diversity_data.csv")
head(d)

d$strategy <- factor(d$strategy, levels = c("CB", "GS-RS"))

## The y-axis is logarithmic, so a ribbon bound of zero or less cannot be drawn.
## mean - sd goes below zero in the last few cycles, so the lower edge is held
## at a small positive floor. Raise or lower 0.5 to taste; on a linear axis
## (scale_y_log10() deleted) you can use plain mean - sd instead.


png("Figure4D_MYMYYY.png", units="in", width=8, height=6, res=600)
d$lo <- pmax(d$mean - d$sd, 0.5)
d$hi <- d$mean + d$sd
p_4D <- ggplot(d, aes(x = year, y = mean, colour = strategy, fill = strategy)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +                     # dots take the line's colour
  scale_y_log10(breaks = c(1, 2, 5, 10, 20, 50, 100)) +
  scale_colour_manual(values = c("CB" = "#2d6d66", "GS-RS" = "#a0522d"),
                      labels = c("Connected Breeding (GS + bridge materials)", "Recurrent genomic selection (closed pool)")) +
  scale_fill_manual(values = c("CB" = "#2d6d66", "GS-RS" = "#a0522d"),
                    labels = c("Connected Breeding (GS + bridge materials)", "Recurrent genomic selection (closed pool)")) +
  labs(x = "", y = "",
       colour = NULL, fill = NULL) +
  theme_classic(base_size = 12) +
  theme(
    ## same coordinate legend as 4C; here it sits lower-left, out of the curves.
    ## c(x, y) are panel fractions -- nudge to move the box.
    legend.position      = c(0.02, 0.10),
    legend.justification = c(0, 0),
    legend.background    = element_rect(fill = "white", colour = NA),
    legend.key           = element_blank(),
    legend.spacing.y     = unit(0, "pt"),
    legend.margin        = margin(2, 4, 2, 2)
  )

p_4D
dev.off()
#ggsave("figures2/Fig_4D.pdf", p_4D, width = 6, height = 4.5)


################################################################################
## FIGURE S2A -- cumulative gain across the donor-inflow ladder
## columns: scenario, label, year, gain, gain_sd, genicPct, genicPct_sd, n
################################################################################

d <- read.csv("figures/Figure_S2_AB_inflow_data.csv")
head(d)

## put the five scenarios in ladder order rather than alphabetical order
#d$label <- factor(d$label,
                #  levels = c("GS-RS (no inflow)", "1 bridge line / cycle",
                            # "2 (as specified)", "4", "8"))

#p_S2A <- ggplot(d, aes(x = year, y = gain, colour = label)) +
 # geom_line(linewidth = 1) +
  #scale_colour_manual(values = c("GS-RS (no inflow)"     = "#a0522d",
                    #             "1 bridge line / cycle" = "#2d6d66",
                        #         "2 (as specified)"      = "#3987e5",
                         #        "4"                     = "#256abf",
                           #      "8"                     = "#104281")) +
  #labs(x = "Year", y = "Cumulative genetic gain (t/ha)",
       #colour = "Bridge lines admitted per cycle") +
 # theme_classic(base_size = 12) +
 # theme(legend.position = "right")
png("FIGURES2_AB_MYYYYY.png", units="in", width=8, height=6, res=600) 
  d$label <- factor(d$label,
                    levels = c("GS-RS (no inflow)", "1 bridge line / cycle",
                               "2 (as specified)", "4", "8"))
  
  pal <- c("GS-RS (no inflow)"     = "#a0522d",
           "1 bridge line / cycle" = "#2d6d66",
           "2 (as specified)"      = "#3987e5",
           "4"                     = "#256abf",
           "8"                     = "#104281")
  
  leg <- c("GS-RS (no inflow)"     = "0  (GS-RS, no inflow)",
           "1 bridge line / cycle" = "1  bridge line",
           "2 (as specified)"      = "2  bridge lines ",
           "4"                     = "4  bridge lines",
           "8"                     = "8  bridge lines")
  
  p_S2A <- ggplot(d, aes(x = year, y = gain, colour = label)) +
    geom_line(linewidth = 1) +
    scale_colour_manual(values = pal,
                        breaks = levels(d$label),
                        labels = leg) +
    labs(x = "Year", y = "Cumulative genetic gain (t/ha)",
         colour = "Bridge lines admitted per cycle") +
    theme_classic(base_size = 12) +
    theme(legend.position = "right")

p_S2A
dev.off()
ggsave("figures2/Fig_S2A.pdf", p_S2A, width = 7.5, height = 4.5)


################################################################################
## FIGURE S2B -- genic variance retained across the same ladder
## same file as S2A; plots genicPct instead of gain
################################################################################

#d <- read.csv("figures/Figure_S2_AB_inflow_data.csv")

#d$label <- factor(d$label,
                #  levels = c("GS-RS (no inflow)", "1 bridge line / cycle",
                            # "2 (as specified)", "4", "8"))

#p_S2B <- ggplot(d, aes(x = year, y = genicPct, colour = label)) +
  #geom_line(linewidth = 1) +
  #scale_y_log10(breaks = c(2, 5, 10, 20, 50, 100)) +
  #scale_colour_manual(values = c("GS-RS (no inflow)"     = "#eb6834",
                                # "1 bridge line / cycle" = "#86b6ef",
                                # "2 (as specified)"      = "#3987e5",
                                # "4"                     = "#256abf",
                                # "8"                     = "#104281")) +
  #labs(x = "Year", y = "Genic variance retained (% of elite base)",
      # colour = "Bridge lines admitted per cycle") +
  #theme_classic(base_size = 12) +
  #theme(legend.position = "right")


d <- read.csv("figures/Figure_S2_AB_inflow_data.csv")


png("FIGURES2_2B_MYYYYY2ADT.png", units="in", width=8, height=6, res=600) 
lvl <- c("GS-RS (no inflow)", "1 bridge line / cycle",
         "2 (as specified)", "4", "8")
d$label <- factor(d$label, levels = lvl)

pal <- c("GS-RS (no inflow)"     = "#a0522d",
         "1 bridge line / cycle" = "#2d6d66",
         "2 (as specified)"      = "#3987e5",
         "4"                     = "#256abf",
         "8"                     = "#104281")

leg <- c("GS-RS (no inflow)"     = "0  (GS-RS, no inflow)",
         "1 bridge line / cycle" = "1  bridge line",
         "2 (as specified)"      = "2  bridge lines (baseline)",
         "4"                     = "4  bridge lines",
         "8"                     = "8  bridge lines")

p_S2B <- ggplot(d, aes(x = year, y = genicPct, colour = label)) +
  geom_line(linewidth = 1) +
  scale_y_log10(breaks = c(2, 5, 10, 20, 50, 100),
                labels = function(x) paste0(x, "%")) +
  scale_colour_manual(values = pal,
                      breaks = lvl,
                      limits = lvl,
                      labels = leg) +
  labs(x = "Year", y = "Genic variance retained (% of elite base)",
       colour = "Bridge lines admitted per cycle") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")

p_S2B
dev.off()

ggsave("figures2/Fig_S2B.pdf", p_S2B, width = 7.5, height = 4.5)


################################################################################
## FIGURE S2C -- rate of gain over the final five cycles
## columns: scenario, label, n, slope, sd, lo, hi
## This is the number that says whether a strategy is still improving or has
## flattened out. Bars are means, whiskers are 95% confidence intervals.
################################################################################


d <- read.csv("figures/Figure_S2_C_late_slope_data.csv")

png("FIGURES2_2C_MYYYYY2ADT.png", units="in", width=8, height=6, res=600) 
pal <- c("GS-RS (no inflow)"     = "#a0522d",
         "1 bridge line / cycle" = "#2d6d66",
         "2 (as specified)"      = "#3987e5",
         "4"                     = "#256abf",
         "8"                     = "#104281")

leg <- c("GS-RS (no inflow)"     = "0  (GS-RS, no inflow)",
         "1 bridge line / cycle" = "1  bridge line",
         "2 (as specified)"      = "2  bridge lines (baseline)",
         "4"                     = "4  bridge lines",
         "8"                     = "8  bridge lines")

## order the bars smallest to largest
d$label <- factor(d$label, levels = d$label[order(d$slope)])

p_S2C <- ggplot(d, aes(x = slope, y = label, fill = label)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y",
                width = 0.2, colour = "grey30") +
  scale_fill_manual(values = pal, breaks = names(pal), guide = "none") +
  scale_y_discrete(labels = leg) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = "Rate of gain, years 45-60 (t/ha/yr)", y = NULL,
       fill = "Bridge lines admitted per cycle") +
  theme_classic(base_size = 12)

p_S2C
dev.off()



d <- read.csv("figures/Figure_S2_C_late_slope_data.csv")
head(d)

## order the bars smallest to largest



d$label <- factor(d$label, levels = d$label[order(d$slope)])
p_S2C <- ggplot(d, aes(x = slope, y = label, fill = label)) +
  geom_col(width = 0.7) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, colour = "grey30") +
  scale_fill_manual(values = c("GS-RS (no inflow)"     = "#eb6834",
                               "1 bridge line / cycle" = "#86b6ef",
                               "2 (as specified)"      = "#3987e5",
                               "4"                     = "#256abf",
                               "8"                     = "#104281"),
                    guide = "none") +
  labs(x = "Rate of gain, years 45-60 (t/ha/yr)", y = NULL) +
  theme_classic(base_size = 12)

p_S2C
dev.off()
ggsave("figures2/Fig_S2C.pdf", p_S2C, width = 7, height = 3.5)


################################################################################
## FIGURE S2D -- year-60 difference for every scenario
## columns: scenario, axis, level, n, delta, sd, lo, hi, p, lab
## Bars right of zero favour Connected Breeding. Whiskers are 95% CIs, so any
## bar whose whisker crosses zero is not statistically distinguishable from it.
################################################################################

d <- read.csv("figures/Figure_S2_D_sensitivity_data.csv")
head(d)
png("FIGURES2_2DD_MYYYYY2ADT.png", units="in", width=12, height=8, res=600) 
d$lab <- factor(d$lab, levels = d$lab[order(d$delta)])

p_S2D <- ggplot(d, aes(x = delta, y = lab, fill = delta > 0)) +
  geom_col(width = 0.7) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, colour = "grey30") +
  geom_vline(xintercept = 0, colour = "black") +
  scale_fill_manual(values = c("TRUE" = "#2a78d6", "FALSE" = "#c0442a"),
                    guide = "none") +
  labs(x = "Year-60 difference vs closed pool (t/ha)", y = NULL) +
  theme_classic(base_size = 11)

p_S2D
dev.off()
ggsave("figures2/Fig_S2D.pdf", p_S2D, width = 9, height = 7)


################################################################################
## NOTES
##
## Change a colour           edit the hex codes in that block's scale_*_manual()
## Change axis limits        add  + coord_cartesian(ylim = c(0, 8))
## Add a title               add  + ggtitle("Cumulative genetic gain")
## Drop the SD ribbon        delete the geom_ribbon() line in the 4C block
## Linear instead of log     delete the scale_y_log10() line
## Bigger text               change base_size in theme_classic()
## Save as PNG               ggsave("figures2/Fig_4C.png", p_4C, width = 6,
##                                  height = 4.5, dpi = 300)
##
## Show only some scenarios, e.g. in Figure S2D:
##     d <- d[d$axis == "resource_matched", ]
##   then rebuild the plot.
##
## To plot something not shown here, the per-cycle output has one row per
## scenario x replicate x cycle:
##     raw <- do.call(rbind, lapply(Sys.glob("out/raw_core*.csv"), read.csv))
##     cb  <- raw[raw$scenario == "CB_base", ]
##     ggplot(cb, aes(year, bridgePoolY)) + stat_summary(fun = mean, geom = "line")
##   useful columns: gain, varG, varGpct, varGenic, varGenicPct, bridgePoolY,
##   donorPoolY, nElig, nBridgeIn, nGeno, nPheno
################################################################################
