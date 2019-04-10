#' @import magrittr
#' @import ggplot2
update_namespace <- function (x) {

}

#' Simulate timecourses
#'
#' @param n number of timecourses to
#' @param Pr_impulse probability of an impulse (versus a sigmoid)
#' @inheritParams fit_timecourse
#'
#' @export
#'
#' @examples
#' simulate_timecourses(n = 20)
simulate_timecourses <- function (n, Pr_impulse = 0.3, timepts = seq(0, 90, by = 5), prior_pars = c("v_sd" = 1.2, "rate_shape" = 2, "rate_scale" = 0.25, "time_shape" = 2, "time_scale" = 25)) {

  model_counts <- rmultinom(1, size = n, prob = c(1-Pr_impulse, Pr_impulse))

  model_parameters <- tibble::tibble(model = c("sigmoid", "impulse"),
                 n = c(model_counts)) %>%
    dplyr::mutate(model_pars = purrr::map2(n, model, simulate_parameters, prior_pars)) %>%
    dplyr::select(-n) %>%
    tidyr::unnest(model_pars)

  model_parameters %>%
    tidyr::nest(-tc_id, -model, .key = "params") %>%
    dplyr::mutate(timecourses = purrr::map2(params, model, fit_timecourse, timepts = timepts))
}

#' Simulate parameters
#'
#' @param n # of parameter sets to simulate
#' @inheritParams validate_priors
#'
#' @export
#'
#' @examples
#' simulate_parameters(n = 20)
simulate_parameters <- function (n, model = "sigmoid",  prior_pars = c("v_sd" = 1.2, "rate_shape" = 2, "rate_scale" = 0.25, "time_shape" = 2, "time_scale" = 25)) {

  stopifnot(class(n) %in% c("numeric", "integer"), n >= 1)

  validate_priors(model, prior_pars)

  if (model %in% c("sigmoid", "impulse")) {
    v_inter <- rnorm(n, mean = 0, sd = prior_pars['v_sd'])
    t_rise <- rgamma(n, shape = prior_pars['time_shape'], scale = prior_pars['time_scale'])
    rate <- rgamma(n, shape = prior_pars['rate_shape'], scale = prior_pars['rate_scale'])
  }

  if (model %in% c("impulse")) {
    v_final <- rnorm(n, mean = 0, sd = prior_pars['v_sd'])
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

#' Validate Priors
#'
#' Check that all required prior parameters are provided and have valid values
#'
#' @param model name of model to simulate
#' @param prior_pars prior parameters for model
validate_priors <- function (model, prior_pars) {

  # validate model type
  valid_models <- c("sigmoid", "impulse")
  stopifnot(class(model) == "character", length(model) == 1)

  if (!(model %in% valid_models)) {
    stop (model, " is an invalid model type; valid model types are: ", paste(valid_models, collapse = ", "))
  }

  # validate priors for model type

}

