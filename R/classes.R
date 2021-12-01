#' @importFrom dplyr %>%
#' @importFrom rlang :=
#' @import tensorflow
#' @import ggplot2
utils::globalVariables(c(
  "sim_fit",
  "timecourses",
  "tc_id",
  "init_id",
  "variable",
  "value",
  "time",
  "fit",
  "model",
  "logLik",
  "sigmoid",
  "impulse",
  "logLik_diff",
  "model_pchisq",
  "model_qchisq",
  "t_rise",
  "t_fall",
  "par",
  "params",
  "parameters",
  "measurements",
  "model_pars",
  "time",
  "v_inter",
  "v_final",
  ".",
  "rate",
  "assymp_type",
  "assymp",
  "t_saturation_start",
  "t_saturation_end",
  "abundance",
  "best_model",
  "fitted_timecourses",
  "loss",
  "v_abs_sum"))

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
#' @param conda_env Conda environment name
#'
#' @export
auto_config_tf <- function (conda_env = "r-tensorflow") {

  checkmate::assertString(conda_env)
  conda_envs <- reticulate::conda_list()

  if ("character" %in% class(conda_envs) ||
    !(conda_env %in% conda_envs$name)) {
      # if no environment exists then create a conda environment w/ TF
    reticulate::conda_create(envname = conda_env)
  }

  # force reticulate to use the specified conda environment
  reticulate::use_condaenv(conda_env, required = TRUE)
  # install TF an TF probability if they don't already exist
  tf_install(conda_env)

  if (!reticulate::py_module_available("tensorflow")) {
    stop ("TensorFlow was not found after installation. This may be because the conda path was not found")
  }
}

tf_install <- function (conda_env) {

  if (!reticulate::py_module_available("tensorflow")) {
    print(paste0("Installing Tensorflow in the ", conda_env, " conda environment"))

    tensorflow::install_tensorflow(
      method = "conda",
      envname = conda_env,
      version = 2.5,
      restart_session = FALSE
    )
  }

  if (!reticulate::py_module_available("tensorflow_probability")) {
    print(paste0("Installing TF probability into the ", conda_env, " conda environment"))

    reticulate::conda_install(
      envname = conda_env,
      packages = "tensorflow-probability==0.13.0",
      pip = TRUE
    )
  }

  return(invisible(0))
}
