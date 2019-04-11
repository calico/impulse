#' @import magrittr
#' @import ggplot2
update_namespace <- function (x) {

}

#' Validate Parameters
#'
#' Check that all required prior parameters are provided and have valid values
#'
#' @param model name of model to simulate
#' @param prior_pars prior parameters for model
validate_priors <- function (model, parameters) {

  # validate model type
  valid_models <- c("sigmoid", "impulse")
  stopifnot(class(model) == "character", length(model) == 1)

  if (!(model %in% valid_models)) {
    stop (model, " is an invalid model type; valid model types are: ", paste(valid_models, collapse = ", "))
  }

  # validate parameters for model type

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

