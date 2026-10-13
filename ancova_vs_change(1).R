# Toy RCT: baseline (month 0) + follow-ups at month 3 and month 4 (uneven spacing).
# Time is categorical, so spacing does not enter the models; it only affects the plot x-axis.
# All contrasts and plotted predictions come from the marginaleffects package (no emmeans).
library(lme4); library(marginaleffects); library(ggplot2); library(tidyr)
set.seed(42)
n <- 150
d <- data.frame(ID = 1:(2*n), arm = factor(rep(c("Control","Treat"), each = n)))
g <- as.numeric(d$arm == "Treat")
u <- rnorm(2*n, 0, 6)
d$y0 <- 50 + u + rnorm(2*n, 0, 6)
d$y1 <- 50 + 0.7*(d$y0-50) + 8*g + rnorm(2*n, 0, 7)   # month 3: true effect 8
d$y2 <- 50 + 0.7*(d$y0-50) + 4*g + rnorm(2*n, 0, 7)   # month 4: true effect 4

long <- pivot_longer(d, c(y0,y1,y2), names_to = "time", values_to = "y")
long$time <- factor(long$time, levels = c("y0","y1","y2"), labels = c("0","3","4"))  # months
long$tx <- factor(ifelse(long$arm=="Treat" & long$time!="0", paste0("T", long$time), "ref"),
                  levels = c("ref","T3","T4"))
long$y0 <- d$y0[match(long$ID, d$ID)]
long <- as.data.frame(long)
fu <- droplevels(subset(long, time != "0"))

m1  <- lmer(y ~ time*arm + (1|ID), long)       # y0+y1+y2, arm main effect kept
m1b <- lmer(y ~ time + tx + (1|ID), long)        # cLDA: baseline difference constrained to 0
m2  <- lmer(y ~ time*arm + y0 + (1|ID), fu)    # multilevel ANCOVA (y0 averaged over its observed distribution)

# ---- Contrasts: Treat - Control at each time (population-level, random effects set to 0) ----
tidy <- function(x, lab) data.frame(model = lab, time = as.character(x$time), estimate = x$estimate,
                  SE = x$std.error, lower.CL = x$conf.low, upper.CL = x$conf.high)
c1  <- avg_comparisons(m1, variables = "arm", by = "time", re.form = NA)
c2  <- avg_comparisons(m2, variables = "arm", by = "time", re.form = NA)
# cLDA: the treatment 'switches on' only at follow-up; contrast tx = T<t> vs ref at time t
c1b <- do.call(rbind, lapply(c("3","4"), function(t) {
  x <- comparisons(m1b, variables = list(tx = c("ref", paste0("T", t))),
                   newdata = datagrid(model = m1b, time = t, tx = "ref"), re.form = NA)
  data.frame(time = t, estimate = x$estimate, std.error = x$std.error, conf.low = x$conf.low, conf.high = x$conf.high) }))
c1b <- rbind(data.frame(time = "0", estimate = 0, std.error = NA, conf.low = NA, conf.high = NA), c1b)
res <- rbind(tidy(c1, "1 y0+y1+y2 (time*arm)"),
             tidy(c1b, "1b cLDA (baseline constrained to 0)"),
             tidy(c2, "2 ANCOVA (y0 covariate)"))
num <- sapply(res, is.numeric); res[num] <- round(res[num], 3)
print(res[order(res$time, res$model), ], row.names = FALSE)

# ---- Predicted means for the plot ----
pp <- function(x, lab) data.frame(model = lab, arm = as.character(x$arm), month = as.numeric(as.character(x$time)),
                                  mean = x$estimate, lo = x$conf.low, hi = x$conf.high)
P1 <- pp(avg_predictions(m1, by = c("time","arm"), re.form = NA), "1 y0+y1+y2 ~ time*arm")
P2 <- pp(avg_predictions(m2, by = c("time","arm"), re.form = NA), "2 y1+y2 ~ time*arm + y0")
# month 0 for ANCOVA panel: observed baseline mean (baseline is a covariate, not modelled)
P2 <- rbind(P2, data.frame(model = P2$model[1], arm = levels(d$arm), month = 0,
                           mean = tapply(d$y0, d$arm, mean), lo = NA, hi = NA))
nd <- rbind(data.frame(time = c("0","3","4"), tx = "ref", arm = "Control"),
            data.frame(time = c("0","3","4"), tx = c("ref","T3","T4"), arm = "Treat"))
nd$time <- factor(nd$time, levels = levels(long$time)); nd$tx <- factor(nd$tx, levels = levels(long$tx))
x1b <- predictions(m1b, newdata = nd, re.form = NA)
P1b <- pp(transform(x1b, time = time), "1b cLDA: baseline constrained to 0")
P <- rbind(P1, P1b, P2)

ggplot(P, aes(month, mean, colour = arm, group = arm)) +
  geom_line(position = position_dodge(.1)) +
  geom_pointrange(aes(ymin = lo, ymax = hi), position = position_dodge(.1), na.rm = TRUE) +
  facet_wrap(~ model, ncol = 1) + scale_x_continuous(breaks = c(0,3,4)) +
  labs(y = "Model-predicted mean outcome", x = "Month",
       title = "Predictions via marginaleffects (months 0, 3, 4)") + theme_bw()
ggsave("/mnt/user-data/outputs/predictions.png", width = 7, height = 8, dpi = 130)
