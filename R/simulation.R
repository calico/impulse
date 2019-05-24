#' Simulate timecourses
#'
#' @param n number of timecourses to
#' @param Pr_impulse probability of an impulse (versus a sigmoid)
#' @inheritParams fit_timecourse
#' @inheritParams estimate_timecourse_params_tf
#' @param measurement_sd gaussian measurement noise to add
#'
#' @export
#'
#' @examples
#' simulate_timecourses(n = 20)
simulate_timecourses <- function (n, Pr_impulse = 0.3, timepts = seq(0, 120, by = 5),
                                  prior_pars = c("v_sd" = 1.2, "rate_shape" = 2, "rate_scale" = 0.25, "time_shape" = 2, "time_scale" = 15),
                                  measurement_sd = 0.2) {

  stopifnot(class(Pr_impulse) == "numeric", length(Pr_impulse) == 1, Pr_impulse >= 0, Pr_impulse <= 1)
  stopifnot(class(measurement_sd) == "numeric", length(measurement_sd) == 1, measurement_sd > 0)

  model_counts <- rmultinom(1, size = n, prob = c(1-Pr_impulse, Pr_impulse))

  model_parameters <- tibble::tibble(model = c("sigmoid", "impulse"),
                                     n = c(model_counts)) %>%
    dplyr::mutate(model_pars = purrr::map2(n, model, simulate_parameters, prior_pars)) %>%
    dplyr::select(-n) %>%
    tidyr::unnest(model_pars)

  model_parameters %>%
    tidyr::nest(-tc_id, -model, .key = "params") %>%
    dplyr::mutate(measurements = purrr::map2(params, model, fit_timecourse, timepts = timepts, fit.label  = "sim_fit"),
                  measurements = purrr::map(measurements, add_noise, measurement_sd = measurement_sd),
                  tc_id = 1:dplyr::n()) %>%
    dplyr::rename(true_model = model)
}

add_noise <- function(measurements, measurement_sd) {
  measurements %>%
    dplyr::mutate(abundance = rnorm(dplyr::n(), mean = sim_fit, sd = measurement_sd))
}

#' Simulate parameters
#'
#' @param n # of parameter sets to simulate
#' @inheritParams estimate_timecourse_params_tf
#'
#' @export
#'
#' @examples
#' simulate_parameters(n = 20)
simulate_parameters <- function (n, model = "sigmoid",  prior_pars = c("v_sd" = 1.2, "rate_shape" = 2, "rate_scale" = 0.25, "time_shape" = 2, "time_scale" = 15)) {

  stopifnot(class(n) %in% c("numeric", "integer"), n >= 1)

  validate_priors(model, prior_pars)

  if (model %in% c("sigmoid", "impulse")) {
    v_inter <- rnorm_just_tails(n, mean = 0, sd = prior_pars['v_sd'], alpha = 0.2)
    t_rise <- rgamma(n, shape = prior_pars['time_shape'], scale = prior_pars['time_scale'])
    rate <- rgamma(n, shape = prior_pars['rate_shape'], scale = prior_pars['rate_scale'])
  }

  if (model %in% c("impulse")) {
    v_change <- rnorm_just_tails(n, mean = 0, sd = prior_pars['v_sd'], alpha = 0.2)
    v_final <- v_inter + v_change
    t_diff <- rgamma(n, shape = prior_pars['time_shape'], scale = prior_pars['time_scale'])
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

rnorm_just_tails <- function(n, mean, sd, alpha = 0.05) {

  stopifnot(length(alpha) == 1, alpha < 0.5)

  runif(n, -1*alpha, 1*alpha) %>%
    {ifelse(. < 0, 1 + ., .)} %>%
    {qnorm(p = ., mean = mean, sd = sd)}
}
