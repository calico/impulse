library(impulse)

set.seed(123)
timecourses <- simulate_timecourses(n = 2, observation_level_noise = 0.5)

measurements <- timecourses %>%
  tidyr::unnest_legacy(measurements) %>%
  # separate by true model
  tidyr::nest_legacy(-true_model, .key = "measurements") %>%
  # fit all models to each timecourse
  tidyr::crossing(tibble::tibble(model = c("sigmoid", "impulse")))

test_that("prior-free minimizing weighted-MSE", {

  least_squares_fits <- measurements %>%
    dplyr::mutate(timecourse_list = purrr::map2(
      measurements,
      model,
      estimate_timecourse_params_tf,
      n_initializations = 200,
      use_prior = FALSE,
      fit_intercept = TRUE
    )) %>%
    dplyr::mutate(best_timecourse_list = purrr::map(timecourse_list, reduce_best_timecourse_params, reduction_type = "loss-min")) %>%
    tidyr::unnest_legacy(purrr::map(best_timecourse_list, ~ tibble::as_tibble(purrr::map(., list))))

  # compare abundances to fitted abundances

  possible_nest_vars <- c("rate", "t_rise", "t_fall", "v_inter", "v_final", "tzero_offset")
  fitted_values <- least_squares_fits %>%
    tidyr::unnest_legacy(parameters) %>%
    dplyr::select(tc_id, model, variable, value) %>%
    tidyr::spread(variable, value) %>%
    tidyr::nest(parameters = dplyr::any_of(possible_nest_vars)) %>%
    dplyr::mutate(fitted_timecourses = purrr::map2(parameters, model,
                                                   fit_timecourse,
                                                   timepts = least_squares_fits$measurements[[1]]$time)) %>%
    tidyr::unnest_legacy(fitted_timecourses)

  abund_vs_fit <- least_squares_fits %>%
    tidyr::unnest_legacy(measurements) %>%
    dplyr::left_join(fitted_values, by = c("tc_id", "model", "time")) %>%
    dplyr::group_by(tc_id, model) %>%
    dplyr::mutate(noise = sqrt(noise^2 / mean(noise^2))) %>%
    dplyr::mutate(ss = noise^-2*(abundance-fit)^2) %>%
    # calculate observed loss (weighted MSE)
    dplyr::summarize(eval_loss = mean(ss), .groups = "drop") %>%
    # compare to logLik that was minimized in TF
    dplyr::left_join(
      least_squares_fits %>%
        tidyr::unnest(loss) %>%
        dplyr::select(tc_id, model, loss),
      by = c("tc_id", "model")
    )

  testthat::expect_equal(
    object = abund_vs_fit$eval_loss,
    expected = as.numeric(abund_vs_fit$loss),
    tolerance = 0.0001
  )

})

test_that("Bayesian impulse minimizes MAP estimate", {

  map_fits <- measurements %>%
    dplyr::mutate(timecourse_list = purrr::map2(
      measurements,
      model,
      estimate_timecourse_params_tf,
      n_initializations = 200,
      use_prior = TRUE,
      fit_intercept = TRUE
    )) %>%
    dplyr::mutate(best_timecourse_list = purrr::map(timecourse_list, reduce_best_timecourse_params, reduction_type = "loss-min")) %>%
    tidyr::unnest_legacy(purrr::map(best_timecourse_list, ~ tibble::as_tibble(purrr::map(., list))))

  # compare abundances to fitted abundances

  possible_nest_vars <- c("rate", "t_rise", "t_fall", "v_inter", "v_final", "tzero_offset")
  fitted_values <- map_fits %>%
    tidyr::unnest_legacy(parameters) %>%
    dplyr::select(tc_id, model, variable, value) %>%
    tidyr::spread(variable, value) %>%
    tidyr::nest(parameters = dplyr::any_of(possible_nest_vars)) %>%
    dplyr::mutate(fitted_timecourses = purrr::map2(parameters, model,
                                                   fit_timecourse,
                                                   timepts = map_fits$measurements[[1]]$time)) %>%
    tidyr::unnest_legacy(fitted_timecourses)

  abund_vs_fit <- map_fits %>%
    tidyr::unnest_legacy(measurements) %>%
    dplyr::left_join(fitted_values, by = c("tc_id", "model", "time")) %>%
    dplyr::mutate(normal_logLik = dnorm(abundance, fit, noise, log = TRUE)) %>%
    dplyr::group_by(tc_id, model) %>%
    # calculate logLik
    dplyr::summarize(eval_logLik = sum(normal_logLik), .groups = "drop") %>%
    # compare to logLik that was minimized in TF
    dplyr::left_join(
      map_fits %>%
        tidyr::unnest(loss) %>%
        dplyr::select(tc_id, model, logLik),
      by = c("tc_id", "model")
    )

  testthat::expect_equal(
    object = abund_vs_fit$eval_logLik,
    expected = as.numeric(abund_vs_fit$logLik),
    tolerance = 0.001
  )
})

test_that("Sigmoid Likelihood is correct", {

  #CONDA_ENV <- "base397"
  #CONDA_PATH <- "/scratch4/alex/python39-for-sean/miniforge3/bin/conda"
  #reticulate::use_condaenv(
  #  CONDA_ENV,
  #  conda = CONDA_PATH,
  #  required=TRUE
  #)

  fitted_kinetics <- impulse::estimate_timecourse_params_tf(
    measurements = example_timecourse,
    prior_pars = c(v_sd = 1.2, rate_shape = 2, rate_scale = 0.25, time_shape = 2, time_scale = 15),
    model = "sigmoid",
    use_prior = TRUE,
    n_initializations = 300L,
    learning_rate = 0.2,
    n_iterations = 1000,
    fit_intercept = TRUE
  )

  # select a consensus parameter set
  consensus_kinetics <- impulse::reduce_best_timecourse_params(fitted_kinetics)

  # calculate logLik

  fitted_timecourses <- fit_timecourse(
    consensus_kinetics$parameters %>%
      dplyr::select(variable, value) %>%
      tidyr::spread(variable, value),
    timepts = unique(example_timecourse$time),
    model = "sigmoid"
  )

  combined_data <- example_timecourse %>%
    dplyr::left_join(fitted_timecourses, by = "time")

  #ggplot(combined_data, aes(x = time, y = abundance)) +
  #  geom_point() +
  #  geom_line(aes(y = fit), color = "red") +
  #  ggtitle("Example timecourse with sigmoid fit") +
  #  theme_minimal()

  timecourse_logLik <- combined_data %>%
    dplyr::mutate(
      normal_logLik = dnorm(abundance, fit, sd = 0.2, log = TRUE)
    ) %>%
    dplyr::summarize(logLik = sum(normal_logLik))

  testthat::expect_equal(
    object = timecourse_logLik$logLik,
    expected = consensus_kinetics$loss$logLik[1],
    tolerance = 0.01
    )

})


