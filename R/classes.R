#' @importFrom dplyr %>%
#' @import ggplot2
update_namespace <- function (x) {

}

#' Validate Parameters
#'
#' Check that all required prior parameters are provided and have valid values
#'
#' @inheritParams estimate_timecourse_params_tf
validate_priors <- function (model, prior_pars) {

  # validate model type
  valid_models <- c("sigmoid", "impulse")
  stopifnot(class(model) == "character", length(model) == 1)

  if (!(model %in% valid_models)) {
    stop (model, " is an invalid model type; valid model types are: ", paste(valid_models, collapse = ", "))
  }

  # validate parameters for model type

}
