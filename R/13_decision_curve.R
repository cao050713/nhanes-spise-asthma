# =============================================================================
# 13_decision_curve.R -- decision curve analysis (Figure 4G)
#
# Net benefit is computed directly from the definition, so the numbers can be
# checked by hand:
#
#     NB(pt) = TP/n - (FP/n) * pt/(1 - pt)
#
# where a participant is treated when the predicted probability exceeds the
# threshold pt. "Treat all" is the same expression with everyone treated and
# gives the familiar prev - (1-prev)*pt/(1-pt); "treat none" is zero.
#
# The two strategies under test are the locked model before recalibration and
# after cross-fitted recalibration, both from R/12.
#
#   Rscript R/13_decision_curve.R
# =============================================================================

source("R/utils.R")
need("dplyr", "ggplot2")

ext <- read_step("external_predictions")

y  <- ext$y
n  <- length(y)
pt <- DCA_THRESHOLDS

# ---------------------------------------------------------------------------
# 1. Net benefit for the three curves
# ---------------------------------------------------------------------------
net_benefit <- function(y, p, pt) {
  vapply(pt, function(t) {
    treated <- p >= t
    tp <- sum(treated & y == 1)
    fp <- sum(treated & y == 0)
    tp / length(y) - (fp / length(y)) * (t / (1 - t))
  }, numeric(1))
}

dca <- data.frame(
  threshold = rep(pt, 3),
  label = rep(c("Original", "Recalibrated", "Treat All"), each = length(pt)),
  net_benefit = c(net_benefit(y, ext$original, pt),
                  net_benefit(y, ext$recalibrated, pt),
                  net_benefit(y, rep(1, n), pt)),
  stringsAsFactors = FALSE
)
dca$label <- factor(dca$label,
                    levels = c("Recalibrated", "Original", "Treat All"))
save_table(dca, "Table4_dca_curve.csv")

# ---------------------------------------------------------------------------
# 2. Net benefit at the prespecified thresholds
# ---------------------------------------------------------------------------
nb <- dca[dca$threshold %in% DCA_REPORTED, ]
nb <- nb[order(nb$label, nb$threshold), ]
print(nb, row.names = FALSE)
save_table(nb, "Table4_net_benefit.csv")

# ---------------------------------------------------------------------------
# 3. Threshold window over which each model beats both defaults
# ---------------------------------------------------------------------------
# "Beats both defaults" means net benefit above the "treat all" curve at the
# same threshold (which is itself above "treat none" whenever prevalence is
# the dominating strategy) and therefore also above zero.
window_one <- function(p, label) {
  nb_model <- net_benefit(y, p, pt)
  nb_all   <- net_benefit(y, rep(1, n), pt)
  ok <- nb_model > pmax(nb_all, 0)
  if (!any(ok)) return(data.frame(label = label, lower = NA_real_, upper = NA_real_))
  data.frame(label = label, lower = min(pt[ok]), upper = max(pt[ok]))
}

window <- rbind(window_one(ext$original, "Original"),
                window_one(ext$recalibrated, "Recalibrated"))
print(window, row.names = FALSE)
save_table(window, "Table4_dca_window.csv")

# ---------------------------------------------------------------------------
# 4. Figure 4G
# ---------------------------------------------------------------------------
p <- ggplot2::ggplot(dca, ggplot2::aes(threshold, net_benefit, colour = label)) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::scale_colour_manual(values = c(Original = "grey55",
                                          Recalibrated = "#219EBC",
                                          `Treat All` = "firebrick")) +
  ggplot2::coord_cartesian(ylim = c(-0.02, max(dca$net_benefit) * 1.05)) +
  ggplot2::labs(x = "Threshold probability", y = "Net benefit",
                title = "Decision curve analysis") +
  theme_paper()
save_fig(p, "Figure4G_dca.png", width = 6, height = 5)

stamp_session("13_decision_curve")
