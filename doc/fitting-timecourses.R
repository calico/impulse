## ----setup---------------------------------------------------------------
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
library(tidyr)
library(purrr)
theme_set(theme_bw() + theme(text = element_text(size = 15)))

library(impulse)

## ----simulate_timecourses, fig.height = 10, fig.width = 10---------------
set.seed(1234)
timecourses <- simulate_timecourses(n = 20)

# plot of timecourses + parameters
timecourses %>%
  tidyr::unnest(measurements) %>%
  tidyr::unite(label, true_model, tc_id) %>%
  ggplot(aes(x = time)) +
  geom_point(aes(y = abundance)) +
  geom_path(aes(y = sim_fit)) +
  facet_wrap(~ label)

## ----estimate_parameters-------------------------------------------------
timecourse_parameters <- timecourses %>%
  unnest(measurements) %>%
  # separate by true model
  nest(-true_model, .key = "measurements") %>%
  # fit all models to each timecourse
  crossing(tibble(model = c("sigmoid", "impulse"))) %>%
  mutate(timecourse_list = map2(measurements, model, estimate_timecourse_params_tf, n_initializations = 50)) %>%
  # reduce each timecourse to best fitting model capturing uncertainty around the best minima
  mutate(best_timecourse_list = map(timecourse_list, reduce_best_timecourse_params, reduction_type = "loss-min")) %>%
  unnest(map(best_timecourse_list, ~ as_tibble(map(., list))))

## ----parameter_comparison, fig.height = 5, fig.width = 8-----------------
estimated_parameters <- timecourse_parameters %>%
  filter(true_model == model) %>%
  unnest(parameters) %>%
  select(tc_id, variable, estimated_value = value)
  
true_parameters <- timecourses %>%
  select(-measurements) %>%
  unnest(params) %>%
  gather(variable, true_value, -true_model, -tc_id) %>%
  filter(!is.na(true_value))

true_parameters %>%
  left_join(estimated_parameters, by = c("tc_id", "variable")) %>%
  ggplot(aes(x = true_value, y = estimated_value)) +
  geom_point() +
  facet_wrap(~ true_model + variable, scale = "free", ncol = 5) +
  geom_abline(intercept = 0, slope = 1, color = "blue")

## ----compare_models------------------------------------------------------
model_comparison <- timecourse_parameters %>%
  unnest(loss) %>%
  impulse_sigmoid_comparison(fdr_cutoff = 0.001) 

timecourses %>%
  select(tc_id, true_model) %>%
  left_join(model_comparison %>%
              select(tc_id, best_model),
            by = "tc_id") %>%
  count(true_model, best_model) %>%
  spread(true_model, n, fill = 0)

## ----timecourse_vis, fig.height = 10, fig.width = 10---------------------
fitted_kinetics <- timecourse_parameters %>%
  unnest(parameters) %>%
  select(tc_id, model, variable, value) %>%
  spread(variable, value) %>%
  nest(-tc_id, .key = "fitted_kinetics")

augmented_timecourses <- timecourses %>%
  left_join(fitted_kinetics, by =  "tc_id") %>%
  left_join(model_comparison, by = "tc_id")

kinetics_plotting(augmented_timecourses, saturation = 0.9, max_time = 150, fit_timepoints = 100)

