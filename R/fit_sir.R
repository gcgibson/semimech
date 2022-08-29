cases <- read.csv("/Users/gcg799/Documents/most_recent_student.csv")# Number of students in bed
cases <- cases[order(cases$week_reporting_date ),]
cases <- cases$cases
cases <- cases[40:60]
plot(cases)
# total count
N <- 5000;
plot(cases)
# times
n_days <- length(cases)
t <- seq(0, n_days, by = 1)
t0 = 0
t <- t[-1]

#initial conditions
i0 <- 1
s0 <- N - i0
r0 <- 0
y0 = c(S = s0, I = i0, R = r0)

# data for Stan
data_sir <- list(n_days = n_days, y0 = y0, t0 = t0, ts = t, N = N, cases = cases)

# number of MCMC steps
niter <- 2000
library(rstan)
model <- stan_model("R/model.stan")

fit_sir_negbin <- sampling(model,
                           data = data_sir,
                           iter = niter,
                           chains = 1,
                           seed = 0)


  smr_pred <- cbind(as.data.frame(summary(
  fit_sir_negbin, pars = "pred_cases", probs = c(0.05, 0.5, 0.95))$summary), t, cases)
colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

ggplot(smr_pred, mapping = aes(x = t)) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = "blue", alpha = 0.35) +
  geom_line(mapping = aes(x = t, y = X50.), color = "blue") +
  geom_point(mapping = aes(y = cases)) +
  labs(x = "Day", y = "Number of students in bed")
