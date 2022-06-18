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
    mutate(best_timecourse_list = map(timecourse_list, reduce_best_timecourse_params, reduction_type = "loss-min")) %>%
    unnest_legacy(map(best_timecourse_list, ~ as_tibble(map(., list))))

  # compare abundances to fitted abundances

  possible_nest_vars <- c("rate", "t_rise", "t_fall", "v_inter", "v_final", "tzero_offset")
  fitted_values <- least_squares_fits %>%
    unnest_legacy(parameters) %>%
    select(tc_id, model, variable, value) %>%
    spread(variable, value) %>%
    tidyr::nest(parameters = any_of(possible_nest_vars)) %>%
    dplyr::mutate(fitted_timecourses = purrr::map2(parameters, model,
                                                   fit_timecourse,
                                                   timepts = map_fits$measurements[[1]]$time)) %>%
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
    mutate(best_timecourse_list = map(timecourse_list, reduce_best_timecourse_params, reduction_type = "loss-min")) %>%
    unnest_legacy(map(best_timecourse_list, ~ as_tibble(map(., list))))

  # compare abundances to fitted abundances

  possible_nest_vars <- c("rate", "t_rise", "t_fall", "v_inter", "v_final", "tzero_offset")
  fitted_values <- map_fits %>%
    unnest_legacy(parameters) %>%
    select(tc_id, model, variable, value) %>%
    spread(variable, value) %>%
    tidyr::nest(parameters = any_of(possible_nest_vars)) %>%
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

  ### TO DO - compare MAP parameters to prior estimates

})


