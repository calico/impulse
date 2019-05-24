## ----setup---------------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
theme_set(theme_bw() + theme(text = element_text(size = 15)))

library(impulse)

## ------------------------------------------------------------------------
sigmoid_impulse_plot <- function(timecourse_parameters) {
  
  fit_timecourse(timecourse_parameters, model = "sigmoid", fit.label = "sigmoid") %>%
    dplyr::left_join(fit_timecourse(timecourse_parameters, model = "impulse", fit.label = "impulse"), by = "time") %>%
    tidyr::gather(eqtn, level, -time) %>%
    ggplot(aes(x = time, y = level, color = eqtn)) +
    geom_path(size = 2)
  
  }
sigmoid_impulse_plot(data_frame(t_rise = 25, rate = 0.25, v_inter = 3, v_final = -3, t_fall = 45))

## ------------------------------------------------------------------------
sigmoid_impulse_plot(data_frame(t_rise = 25, rate = -0.25, v_inter = 3, v_final = -3, t_fall = 45)) +
  ggtitle("Negative rate pathology")

## ----fig.height = 8, fig.width = 8---------------------------------------
list(sigmoid_impulse_plot(data_frame(t_rise = 25, rate = 0.1, v_inter = 3, v_final = -3, t_fall = 45)) +
       ggtitle("rate = 0.1"),
     sigmoid_impulse_plot(data_frame(t_rise = 25, rate = 0.5, v_inter = 3, v_final = -3, t_fall = 45)) +
       ggtitle("rate = 0.5"),
     sigmoid_impulse_plot(data_frame(t_rise = 25, rate = 1, v_inter = 3, v_final = -3, t_fall = 45)) +
       ggtitle("rate = 1"),
     sigmoid_impulse_plot(data_frame(t_rise = 25, rate = 3, v_inter = 3, v_final = -3, t_fall = 45)) +
       ggtitle("rate = 0.3")) %>%
  {do.call(grid.arrange, .)}

## ----fig.height = 10, fig.width = 10-------------------------------------
betas <- seq(0, 2, by = 0.01)
shape_scale_prod <- 0.5
shape_scale_ratio <- 10^(seq(0, 2, length.out = 20))
gamma_densities <- lapply(shape_scale_ratio, function(a_shape_scale_ratio) {
  shape = sqrt(shape_scale_prod) * sqrt(a_shape_scale_ratio)
  scale = sqrt(shape_scale_prod) / sqrt(a_shape_scale_ratio)
  data_frame(beta = betas, 
             shape_scale = paste(round(shape, 3), round(scale, 3), sep = "-"),
             density = dgamma(betas, shape = shape, scale = scale))
}) %>%
  dplyr::bind_rows()

ggplot(gamma_densities, aes(x = beta, y = density)) +
  facet_wrap(~ shape_scale) +
  geom_path() +
  ggtitle("shape & scale, where shape * scale = 0.5")

## ----fig.height = 10, fig.width = 6--------------------------------------

gamma_quantiles <- qgamma(c(0.025, 0.5, 0.95), shape = 2, scale = 0.25)
list(sigmoid_impulse_plot(data_frame(t_rise = 25, rate = 0.06, v_inter = 3, v_final = -3, t_fall = 45)) +
       ggtitle("2.5% rate quantile"),
     sigmoid_impulse_plot(data_frame(t_rise = 25, rate = 0.42, v_inter = 3, v_final = -3, t_fall = 45)) +
       ggtitle("50% rate quantile"),
     sigmoid_impulse_plot(data_frame(t_rise = 25, rate = 1.19, v_inter = 3, v_final = -3, t_fall = 45)) +
       ggtitle("97.5% rate quantile")) %>%
       {do.call(grid.arrange, .)}

## ------------------------------------------------------------------------
sigmoid_impulse_plot(data_frame(t_rise = -25, rate = 0.25, v_inter = 3, v_final = -3, t_fall = 90)) +
  ggtitle(expression(t[rise] < 0 ~ " pathology"))

## ------------------------------------------------------------------------
list(sigmoid_impulse_plot(data_frame(t_rise = 25, rate = 0.5, v_inter = 3, v_final = -3, t_fall = 45)) +
       ggtitle(expression(t[rise] ~ "= 25," ~ t[fall] ~ "= 45")),
     sigmoid_impulse_plot(data_frame(t_rise = 45, rate = 0.5, v_inter = 3, v_final = -3, t_fall = 25)) +
       ggtitle(expression(t[rise] ~ "= 45," ~ t[fall] ~ "= 45 pathology")),
     sigmoid_impulse_plot(data_frame(t_rise = 45, rate = 0.5, v_inter = 3, v_final = -3, t_fall = 45)) +
       ggtitle(expression(t[rise] ~ "= 45" ~ t[fall] ~ "= 45 pathology"))) %>%
       {do.call(grid.arrange, .)}

## ------------------------------------------------------------------------
qplot(y = dgamma(seq(0, 200, by = 0.1), shape = 2, scale = 25), x = seq(0, 200, by = 0.1)) +
  scale_x_continuous("time") + scale_y_continuous(expression("Pr(t |" ~ Gamma ~ ")"))

## ----fig.height = 10, fig.width = 10-------------------------------------
gamma_quantiles <- qgamma(c(0.025, 0.25, 0.5, 0.75, 0.95), shape = 2, scale = 25)
list(sigmoid_impulse_plot(data_frame(t_rise = gamma_quantiles[1], rate = 0.5, v_inter = 3, v_final = -3, t_fall = gamma_quantiles[1]*2)) +
       ggtitle("2.5% t quantile"),
     sigmoid_impulse_plot(data_frame(t_rise = gamma_quantiles[2], rate = 0.5, v_inter = 3, v_final = -3, t_fall = gamma_quantiles[2]*2)) +
       ggtitle("25% t quantile"),
     sigmoid_impulse_plot(data_frame(t_rise = gamma_quantiles[3], rate = 0.5, v_inter = 3, v_final = -3, t_fall = gamma_quantiles[3]*2)) +
       ggtitle("50% t quantile"),
     sigmoid_impulse_plot(data_frame(t_rise = gamma_quantiles[4], rate = 0.5, v_inter = 3, v_final = -3, t_fall = gamma_quantiles[4]*2)) +
       ggtitle("75% t quantile")) %>%
       {do.call(grid.arrange, .)}

