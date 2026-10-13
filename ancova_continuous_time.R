# Same toy RCT (months 0, 3, 4; true effect 8 at month 3, 4 at month 4), but time is CONTINUOUS (linear per arm).
library(lme4); library(marginaleffects); library(ggplot2); library(tidyr)
set.seed(42)
n <- 150
d <- data.frame(ID = 1:(2*n), arm = factor(rep(c("Control","Treat"), each = n)))
g <- as.numeric(d$arm == "Treat")
u <- rnorm(2*n, 0, 6)
d$y0 <- 50 + u + rnorm(2*n, 0, 6)
d$y1 <- 50 + 0.7*(d$y0-50) + 8*g + rnorm(2*n, 0, 7)
d$y2 <- 50 + 0.7*(d$y0-50) + 4*g + rnorm(2*n, 0, 7)

long <- as.data.frame(pivot_longer(d, c(y0,y1,y2), names_to = "v", values_to = "y"))
long$month <- c(y0 = 0, y1 = 3, y2 = 4)[long$v]
long$y0 <- d$y0[match(long$ID, d$ID)]
fu <- subset(long, month > 0)

# Categorical-time reference: saturated means (what the earlier script estimated)
mcat <- lmer(y ~ factor(month) + I(arm=="Treat" & month==3) + I(arm=="Treat" & month==4) + (1|ID), long)

c1  <- lmer(y ~ month*arm + (1|ID), long)    # lines per arm, free baseline difference
c1b <- lmer(y ~ month + month:arm + (1|ID), long)  # cLDA: common intercept at month 0, one slope per arm
c2  <- lmer(y ~ month*arm + y0 + (1|ID), fu) # ANCOVA, follow-ups only (2 points per arm => line is saturated)

tidy <- function(x, lab) data.frame(model = lab, month = x$month, estimate = x$estimate,
                  SE = x$std.error, lower.CL = x$conf.low, upper.CL = x$conf.high)
con <- function(m, lab, months = c(0,3,4)) tidy(avg_comparisons(m, variables = "arm", by = "month",
                  newdata = subset(if (identical(m, c2)) fu else long, month %in% months), re.form = NA), lab)
res <- rbind(con(c1, "C1  y0+y1+y2 ~ month*arm"),
             con(c1b, "C1b cLDA, linear (common intercept)"),
             con(c2, "C2  ANCOVA, linear time"))
res$truth <- c(0,8,4)[match(res$month, c(0,3,4))]
num <- sapply(res, is.numeric); res[num] <- round(res[num], 3)
print(res[order(res$month, res$model), ], row.names = FALSE)
cat("\nLinear-time slope difference (Treat - Control, per month):\n")
print(round(coef(summary(c1))["month:armTreat", ], 3)); print(round(coef(summary(c1b))[grep("armTreat", rownames(coef(summary(c1b)))), ], 3))

# ---- Plot: fitted straight lines + observed means ----
grid <- function(m, months) datagrid(model = m, month = months, arm = c("Control","Treat"))
pl <- function(m, lab, months) { x <- predictions(m, newdata = grid(m, months), re.form = NA)
  data.frame(model = lab, arm = as.character(x$arm), month = x$month, mean = x$estimate, lo = x$conf.low, hi = x$conf.high) }
P <- rbind(pl(c1,  "C1  y0+y1+y2 ~ month*arm",           seq(0,4,.25)),
           pl(c1b, "C1b cLDA, linear (common intercept)", seq(0,4,.25)),
           pl(c2,  "C2  ANCOVA, linear time",             seq(3,4,.25)))
obs <- aggregate(y ~ arm + month, long, mean)
O <- do.call(rbind, lapply(unique(P$model), function(m) cbind(model = m, obs)))
ggplot(P, aes(month, mean, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, colour = NA) + geom_line() +
  geom_point(data = O, aes(y = y), shape = 21, fill = "white", size = 2.5, stroke = 1) +
  facet_wrap(~ model, ncol = 1) + scale_x_continuous(breaks = c(0,3,4)) +
  labs(y = "Model-predicted mean outcome", x = "Month",
       title = "Continuous time: fitted lines (open circles = observed means)") + theme_bw()
ggsave("/mnt/user-data/outputs/predictions_continuous.png", width = 7, height = 8, dpi = 130)
