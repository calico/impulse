#' @importFrom dplyr %>%
#' @importFrom rlang :=
#' @import tensorflow
#' @import ggplot2
utils::globalVariables(c("sim_fit", "timecourses", "tc_id", "init_id",
                         "variable", "value", "time", "fit", "model", "logLik",
                         "sigmoid", "impulse", "logLik_diff", "model_pchisq",
                         "model_qchisq", "t_rise", "t_fall", "par", "params",
                         "parameters", "measurements", "model_pars", "time",
                         "v_inter", "v_final", ".", "rate", "assymp_type",
                         "assymp", "t_saturation_start", "t_saturation_end",
                         "abundance", "best_model", "fitted_timecourses",
                         "loss", "v_abs_sum"))

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
    stop (model, " is an invalid model type; valid model types are: ",
          paste(valid_models, collapse = ", "))
  }

  # validate parameters for model type

  stopifnot(class(prior_pars) == "numeric")
  required_pars <- c("v_sd", "rate_shape", "rate_scale", "time_shape",
                     "time_scale")
  provided_pars <- names(prior_pars)
  missing_required_pars <- setdiff(required_pars, provided_pars)

  if (length(missing_required_pars) != 0) {
    stop ("missing ", length(missing_required_pars),
          " required parameters in \"prior_pars\": ",
          paste(missing_required_pars, collapse = ", "))
  }

  stopifnot(prior_pars["v_sd"] > 0)
  stopifnot(prior_pars["rate_shape"] > 0)
  stopifnot(prior_pars["rate_scale"] > 0)
  stopifnot(prior_pars["time_shape"] > 0)
  stopifnot(prior_pars["time_scale"] > 0)
}

#' Auto Configure TensforFlow
#'
#' Load and install the r-tensorflow conda environment (python / tensorflow
#' can be setup in other ways with reticulate).
#'
#' @export
auto_config_tf <- function () {

  conda_envs <- reticulate::conda_list()

  if ("character" %in% class(conda_envs) ||
      !("r-tensorflow" %in% conda_envs$name)) {
    tensorflow::install_tensorflow(method = "conda",
                                   envname = "r-tensorflow",
                                   version = 2.4,
                                   extra_packages = "tensorflow-probability")
  } else {
    reticulate::use_condaenv("r-tensorflow", required = TRUE)

    if (!reticulate::py_module_available("tensorflow")) {
      tensorflow::install_tensorflow(method = "conda",
                                     envname = "r-tensorflow",
                                     version = 2.4,
                                     extra_packages = "tensorflow-probability")
    }
  }

  if (!reticulate::py_module_available("tensorflow")) {
    stop ("TensorFlow was not foound after installation")
  }

  tf_v1_compatibility()
}

tf_v1_compatibility <- function () {

  library(tensorflow)
  #tensorflow::use_compat(version = "v1")
  tf$compat$v1$disable_eager_execution()

}
