#' Simulate timecourses
#'
#' @param n number of timecourses to
#' @param Pr_impulse probability of an impulse (versus a sigmoid)
#' @inheritParams fit_timecourse
#' @inheritParams estimate_timecourse_params_tf
#' @param measurement_sd gaussian measurement noise to add
#' @param observation_level_noise scale samples by
#'   exp(rnorm(1)*\code{observation_level_noise}) to introduce samples with
#'   varying noise levels. Set a zero by default for a homoschedastic model.
#'
#' @export
#'
#' @examples
#' simulate_timecourses(n = 20)
simulate_timecourses <- function (
    n,
    Pr_impulse = 0.5,
    timepts = seq(0, 120, by = 5),
    prior_pars = c("v_sd" = 1.2, "rate_shape" = 2, "rate_scale" = 0.25,
                   "time_shape" = 1, "time_scale" = 30),
    measurement_sd = 0.2,
    observation_level_noise = 0
    ) {

  checkmate::assertNumber(Pr_impulse, lower = 0, upper = 1)
  checkmate::assertNumber(measurement_sd, lower = 0)
  checkmate::assertNumber(observation_level_noise, lower = 0)

  model_counts <- stats::rmultinom(1,
                                   size = n,
                                   prob = c(1 - Pr_impulse, Pr_impulse))

  model_parameters <- tibble::tibble(model = c("sigmoid", "impulse"),
                                     n = c(model_counts)) %>%
    dplyr::filter(n > 0) %>%
    dplyr::mutate(model_pars = purrr::map2(n, model,
                                           simulate_parameters,
                                           prior_pars)) %>%
    dplyr::select(-n) %>%
    tidyr::unnest_legacy(model_pars)

  timecourses <- model_parameters %>%
    tidyr::nest_legacy(-tc_id, -model, .key = "params") %>%
    dplyr::mutate(
      measurements = purrr::map2(params, model,
                                 fit_timecourse,
                                 timepts = timepts,
                                 fit.label  = "sim_fit"),
      measurements = purrr::map(measurements,
                                add_noise,
                                measurement_sd = measurement_sd,
                                observation_level_noise = observation_level_noise),
      tc_id = 1:dplyr::n()) %>%
    dplyr::rename(true_model = model)

  return(timecourses)
}

add_noise <- function(measurements, measurement_sd, observation_level_noise) {
  measurements %>%
    dplyr::mutate(
      noise = exp(rnorm(dplyr::n())*observation_level_noise)*measurement_sd,
      abundance = stats::rnorm(dplyr::n(),
                               mean = sim_fit,
                               sd = noise)
      )
}

#' Simulate parameters
#'
#' @param n # of parameter sets to simulate
#' @inheritParams estimate_timecourse_params_tf
#' @param alpha tail fraction to sample from for asymptote changes
#'
#' @export
#'
#' @examples
#' simulate_parameters(n = 20)
simulate_parameters <-
  function (n, model = "sigmoid",
            prior_pars = c("v_sd" = 1.2, "rate_shape" = 2, "rate_scale" = 0.25,
                           "time_shape" = 5, "time_scale" = 5),
            alpha = 0.01) {

  stopifnot(class(n) %in% c("numeric", "integer"),
            n >= 0)

  validate_priors(model, prior_pars)

  if (n == 0) {
    return (NULL)
  }

  if (model %in% c("sigmoid", "impulse")) {
    v_inter <- rnorm_just_tails(n,
                                mean_val = 0,
                                sd_val = prior_pars["v_sd"],
                                alpha = alpha)
    t_rise <- stats::rgamma(n,
                            shape = prior_pars["time_shape"],
                            scale = prior_pars["time_scale"])
    rate <- stats::rgamma(n,
                          shape = prior_pars["rate_shape"],
                          scale = prior_pars["rate_scale"])
  }

  if (model %in% c("impulse")) {
    v_change <- rnorm_just_tails(n,
                                 mean_val = 0,
                                 sd_val = prior_pars["v_sd"],
                                 alpha = alpha)
    v_final <- v_inter + v_change
    t_diff <- stats::rgamma(n,
                            shape = prior_pars["time_shape"],
                            scale = prior_pars["time_scale"])
    t_fall <- t_rise + t_diff
  }

  if (model == "sigmoid") {
    tibble::tibble(tc_id = seq(n),
                   v_inter = v_inter,
                   t_rise = t_rise,
                   rate = rate)
  } else if (model == "impulse") {
    tibble::tibble(tc_id = seq(n),
                   v_inter = v_inter,
                   t_rise = t_rise,
                   v_final = v_final,
                   t_fall = t_fall,
                   rate = rate)
  } else {
    stop ("required parameters for \"", model, "\" not defined")
  }
}

rnorm_just_tails <- function(n, mean_val, sd_val, alpha = 0.01) {

  stopifnot(length(alpha) == 1, alpha < 0.5)

  stats::runif(n, -1 * alpha, 1 * alpha) %>%
    {ifelse(. < 0, 1 + ., .)} %>%
    {stats::qnorm(p = ., mean = mean_val, sd = sd_val)}
}
