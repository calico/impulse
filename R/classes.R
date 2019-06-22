#' @importFrom dplyr %>%
#' @import tensorflow
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

  stopifnot(class(prior_pars) == "numeric")
  required_pars <- c("v_sd", "rate_shape", "rate_scale", "time_shape", "time_scale")
  provided_pars <- names(prior_pars)
  missing_required_pars <- setdiff(required_pars, provided_pars)

  if (length(missing_required_pars) != 0) {
    stop ("missing ", length(missing_required_pars), " required parameters in \"prior_pars\": ", paste(missing_required_pars, collapse = ", "))
  }

  stopifnot(prior_pars['v_sd'] > 0)
  stopifnot(prior_pars['rate_shape'] > 0)
  stopifnot(prior_pars['rate_scale'] > 0)
  stopifnot(prior_pars['time_shape'] > 0)
  stopifnot(prior_pars['time_scale'] > 0)
}

#' Auto Configurate TensforFlow
#'
#' Load and install the r-tensorflow conda environment (python / tensorflow can be setup in other ways with reticulate).
#'
#' @export
auto_config_tf <- function () {

  conda_envs <- reticulate::conda_list()

  if ("character" %in% class(conda_envs) || !("r-tensorflow" %in% conda_envs$name)) {
    tensorflow::install_tensorflow()
  } else {
    reticulate::use_condaenv("r-tensorflow")

    if (!reticulate::py_module_available("tensorflow")) {
      tensorflow::install_tensorflow()
    }
  }
}
