# Simulation of the repressilator


library(deSolve)
library(ggplot2)
library(tidyr)

Elowitz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dm1 <- alpha0 + alpha / (1 + p3^n) - m1
    dm2 <- alpha0 + alpha / (1 + p1^n) - m2
    dm3 <- alpha0 + alpha / (1 + p2^n) - m3
    dp1 <- -beta * (p1 - m1)
    dp2 <- -beta * (p2 - m2)
    dp3 <- -beta * (p3 - m3)
    list(c(dm1, dm2, dm3, dp1, dp2, dp3))
  })
}

parameters <- c(n = 2.0, alpha0 = 0.2, alpha = 200.0, beta = 3.0)
state      <- c(m1 = 100, m2 = 80, m3 = 50, p1 = 10, p2 = 10, p3 = 10)
times      <- seq(0, 100, by = 0.1)

out <- ode(y = state, times = times, func = Elowitz, parms = parameters)
df <- data.frame(
  time = out[,"time"],
  lacl = out[,"p1"],
  tetR = out[,"p2"],
  cl   = out[,"p3"]
)

plotElowitz <-
  df %>%
  gather("protein", "value", 2:4) %>%
  ggplot(aes(x=time, y=value, group=protein, colour=protein)) +
  geom_line() +
  theme_bw() +
  ggtitle("Repressilator") +
  xlab("Time (minutes)") +
  ylab("Proteins per cell")

plotElowitz
