#' Fit Timecourses
#'
#' @param measurements a tibble containing:
#' \describe{
#'   \item{tc_id}{a unique indicator for each timecourse"}
#'   \item{time}{a numeric predictor variable}
#'   \item{abundance}{a numeric response variable}
#'   \item{noise (optional)}{if provided, a numeric variable specifying the
#'   noise estimate for each observation. This will be used to calculate the
#'   Gaussian likelihood if priors are used, or to calculate a weighted
#'   sum-of-squares estimate if priors are not used}
#'   }
#' @param model model to fit:
#' \describe{
#'   \item{sigmoid}{one sigmoidal response}
#'   \item{impulse}{two sigmoidal responses}
#'   }
#' @param n_initializations Number of initializations to use for each timecourse.
#' @param max_n_initializations Maximum number of initializations that can be used
#' @param use_prior If FALSE, fit least squares. If TRUE, add priors for a MAP estimate.
#' @param prior_pars Named numeric vector of parameters to use for priors (if use_prior is TRUE)
#' \describe{
#'   \item{v_sd}{Gaussian SD of assymptotes (v_inter and v_final)}
#'   \item{rate_shape}{Shape of rate Gamma}
#'   \item{rate_scale}{Scale of rate Gamma}
#'   \item{time_shape}{Shape of t_rise and t_fall - t_rise Gamma}
#'   \item{time_scale}{Scale of t_rise and t_fall - t_rise Gamma}
#'   }
#' @param fit_intercept If TRUE, the intercept will be fit, if FALSE, the intercept will be constrainted to zero
#' @param learning_rate learning rate for the Adams optimizer
#' @param verbose if TRUE then print additional information
#'
#' @return a timecourse list:
#' \describe{
#'   \item{invalid_timecourse_fits}{tibble of parameter initializations for initializations
#'   that went to NaN for debugging}
#'   \item{loss}{tibble of losses for each tc_id and inititalization (init_id)}
#'   \item{parameters}{tibble of final parameters for each tc_id and initialization (init_id)}
#'   }
#'
#' @examples
#'
#' library(dplyr)
#' library(impulse)
#' auto_config_tf()
#'
#' set.seed(123)
#' timecourses <- simulate_timecourses(n = 2, observation_level_noise = 0.5)
#'
#' z <- timecourses %>%
#'   tidyr::unnest_legacy(measurements) %>%
#'   # separate by true model
#'   tidyr::nest_legacy(-true_model, .key = "measurements") %>%
#'   # fit all models to each timecourse
#'   tidyr::crossing(tibble::tibble(model = c("sigmoid", "impulse"))) %>%
#'   dplyr::mutate(timecourse_params = purrr::map2(
#'     measurements,
#'     model,
#'     estimate_timecourse_params_tf,
#'     n_initializations = 20,
#'     fit_intercept = TRUE
#'     ))
#'
#' @export
estimate_timecourse_params_tf <- function(
  measurements,
  model = "sigmoid",
  n_initializations = 100,
  max_n_initializations = 1000,
  use_prior = TRUE,
  prior_pars = c(
    "v_sd" = 1.2,
    "rate_shape" = 2,
    "rate_scale" = 0.25,
    "time_shape" = 2,
    "time_scale" = 15),
  fit_intercept = FALSE,
  learning_rate = 0.1,
  verbose = FALSE
  ) {

  if (!requireNamespace("tensorflow", quietly = TRUE)) {
    stop('The "tensorflow" package must be installed to use this function',
         call. = FALSE)
  } else {
    # use the v1 tensorflow api
    tf$compat$v1$disable_eager_execution()
    # import probability module
    tfp <- reticulate::import("tensorflow_probability")
  }

  stopifnot("data.frame" %in% class(measurements))
  required_vars <- c("tc_id", "time", "abundance")
  missing_vars <- setdiff(required_vars, colnames(measurements))
  if (length(missing_vars) != 0) {
    stop ("required variables are missing from \"measurements\": ",
          paste(missing_vars, collapse = ", "))
  }

  checkmate::assertChoice(model, c("sigmoid", "impulse"))
  checkmate::assertNumber(n_initializations, lower = 10)
  n_initializations <- as.integer(n_initializations)
  checkmate::assertNumber(max_n_initializations, lower = n_initializations)
  checkmate::assertLogical(use_prior, len = 1)
  checkmate::assertLogical(fit_intercept, len = 1)
  checkmate::assertNumber(learning_rate, lower = 0)
  checkmate::assertLogical(verbose, len = 1)

  # test parameters supplied for parameter initialization / priors

  if (use_prior) {
    stopifnot(length(prior_pars) == length(unique(names(prior_pars))))

    missing_pars <- setdiff(c("v_sd", "rate_shape", "rate_scale", "time_shape",
                              "time_scale"), names(prior_pars))

    if (length(missing_pars) != 0) {
      stop("\"use_prior\" is TRUE, but ", length(missing_pars),
           " required parameters are missing - supply ",
           paste(missing_pars, collapse = ", "), " with \"prior_pars\"")
    }
  } else {
    prior_pars <- c(
      "v_sd" = stats::sd(timecourses$log2_fc),
      "t_max" = max(timecourses$time)
      )
  }

  all_timecourse_fits <- list()
  entry_number <- 0
  current_n_initializations <- n_initializations

  for (a_tc_id in unique(measurements$tc_id)) {
    entry_number <- entry_number + 1

    if (verbose) {
      print(glue::glue("timecourse {a_tc_id} is now running"))
    }

    one_timecourse <- measurements %>%
      dplyr::filter(tc_id == a_tc_id)

    continue <- TRUE
    while (continue) {

      timecourse_fit <- estimate_one_timecourse_params_tf(
        one_timecourse = one_timecourse,
        a_tc_id = a_tc_id,
        use_prior = use_prior,
        n_initializations = current_n_initializations,
        model = model,
        prior_pars = prior_pars,
        fit_intercept = fit_intercept,
        learning_rate = learning_rate,
        verbose = verbose
      )

      if (is.null(timecourse_fit)) {
        current_n_initializations <- current_n_initializations * 2

        if (current_n_initializations > max_n_initializations) {
          warning (glue::glue(
            "The number of initialization has reached the limit set by max_n_initializations: {max_n_initializations}
            for tc_id {a_tc_id}. Fitting will continue with this number of intializations"
            ))
        } else {
          message (glue::glue(
            "Too few valid initializations for tc_id: {a_tc_id}, increasing initializations to {current_n_initializations}"
          ))
        }
      } else {
        continue <- FALSE
      }
    }

    all_timecourse_fits[[entry_number]] <- timecourse_fit
  }

  all_timecourse_fits %>%
    purrr::transpose() %>%
    purrr::map(dplyr::bind_rows)
}

#' Estimate One Timecourses Parameters
#'
#' @param one_timecourse dynamics of a single timecourse including
#'   time and abundance variables
#' @param a_tc_id timecourse id for the timecourse being fit
#' @inheritParams estimate_timecourse_params_tf
estimate_one_timecourse_params_tf <- function (
  one_timecourse,
  a_tc_id,
  use_prior,
  n_initializations,
  model,
  prior_pars = NULL,
  fit_intercept = FALSE,
  learning_rate = 0.1,
  verbose = FALSE
  ) {

  checkmate::assertDataFrame(one_timecourse)
  stopifnot(c("time", "abundance") %in% colnames(one_timecourse))
  checkmate::assertLogical(use_prior, len = 1)
  checkmate::assertNumber(n_initializations, lower = 10)
  checkmate::assertNumber(learning_rate, lower = 0)
  checkmate::assertLogical(verbose, len = 1)
  checkmate::assertLogical(fit_intercept, len = 1)

  tfp <- reticulate::import("tensorflow_probability")

  # Setup initialization

  if (model %in% c("sigmoid", "impulse")) {
    # for parameters shared by sigmoid and impulse

    sigmoid_params <- init_sigmoid_parameters(use_prior, n_initializations, prior_pars)
    v_inter <- sigmoid_params$v_inter
    t_rise <- sigmoid_params$t_rise
    rate <- sigmoid_params$rate
    parameters <- c("v_inter", "t_rise", "rate")
  }

  if (model == "impulse") {
    # setup impulse specific parameters

    impulse_params <- init_impulse_parameters(use_prior, n_initializations, prior_pars, t_rise)
    v_final <- impulse_params$v_final
    t_diff <- impulse_params$t_diff
    t_fall <- impulse_params$t_fall
    parameters <- c(parameters, "v_final", "t_fall")
  }

  # Setup model

  # data
  timepts <- tf$compat$v1$placeholder(
    tf$float32,
    shape(NULL, n_initializations),
    name = "time"
  )
  expression <- tf$compat$v1$placeholder(
    tf$float32,
    shape(NULL, n_initializations),
    name = "measured_expression"
  )
  noise <- tf$compat$v1$placeholder(
    tf$float32,
    shape(NULL, n_initializations),
    name = "noise"
  )

  # either fix offset to zero or keep it as a free parameter
  if (fit_intercept) {
    tzero_offset <- tf$Variable(
      tf$random$normal(shape(n_initializations),
                       mean = 0,
                       stddev = 0),
      name = "tzero_offset")
    parameters <- c(parameters, "tzero_offset")
  } else {
    tzero_offset <- 0
  }

  if (model == "sigmoid") {
    # "v_inter*(1/(1 + exp(-1*rate*(time - t_rise))))"
    # special case of the impulse where h1 = h2

    rise_exp <- tf$exp(tf$multiply(-1, rate) * tf$subtract(timepts, t_rise),
                       name = "exponentiation")

    rise_sigmoid <- tf$multiply(v_inter, tf$divide(1, 1 + rise_exp),
                                name = "expression_sigmoid")

    fit_expression <- tf$add(rise_sigmoid, tzero_offset, name = "expression_fitted_values")

  } else if (model == "impulse") {
    # "(1/(1 + exp(-1*rate*(time - t_rise)))) * (v_final + (v_inter -
    # v_final)*(1/(1 + exp(rate*(time - t_fall)))))"

    rise_exp <- tf$exp(tf$multiply(-1, rate) * tf$subtract(timepts, t_rise),
                       name = "rise_exponentiation")
    rise_sigmoid <- tf$divide(1, 1 + rise_exp, name = "rise_sigmoid")

    fall_exp <- tf$exp(tf$multiply(rate, tf$subtract(timepts, t_fall)),
                       name = "fall_exponentiation")
    fall_sigmoid <- fall_sigmoid <- tf$divide(1, 1 + fall_exp,
                                              name = "fall_sigmoid")
    offset_fall_sigmoid <- v_final +
      (tf$subtract(v_inter, v_final) * fall_sigmoid)

    impulse_dynamics <- tf$multiply(rise_sigmoid, offset_fall_sigmoid,
                                    name = "impulse_dynamics")

    fit_expression <- tf$add(impulse_dynamics, tzero_offset, name = "expression_fitted_values")

  } else {
    stop('model must be either "sigmoid" or "impulse"')
  }

  # Setup priors

  if (use_prior) {
    v_prior <- tfp$distributions$Normal(loc = 0,
                                        scale = prior_pars["v_sd"])
    rate_prior <- tfp$distributions$Gamma(
      concentration = prior_pars["rate_shape"],
      rate = 1 / prior_pars["rate_scale"])
    time_prior <- tfp$distributions$Gamma(
      concentration = prior_pars["time_shape"],
      rate = 1 / prior_pars["time_scale"])

    if (model == "sigmoid") {
      model_log_pr <- tf$add(tf$add(v_prior$log_prob(v_inter),
                                    rate_prior$log_prob(rate)),
                             time_prior$log_prob(t_rise))
    } else if (model == "impulse") {
      model_log_pr <- tf$add(tf$add(tf$add(tf$add(v_prior$log_prob(v_inter),
                                                  rate_prior$log_prob(rate)),
                                           time_prior$log_prob(t_rise)),
                                    v_prior$log_prob(v_final)),
                             time_prior$log_prob(t_diff))
    } else {
      stop('model must be either "sigmoid" or "impulse"')
    }
  }

  # general model formatting
  optimizer <- tf$compat$v1$train$AdamOptimizer(learning_rate)

  if (use_prior) {
    # minimize normal likelihood with priors
    norm_target <- tfp$distributions$Normal(loc = expression,
                                            scale = noise)
    normal_logLik <- tf$reduce_sum(norm_target$log_prob(fit_expression),
                                   axis = 0L, name = "normal_logLik")
    logPr <- tf$subtract(0, tf$add(normal_logLik, model_log_pr))

    # minimize negative logLik + logPrior (max logLik)
    loss <- tf$reduce_sum(logPr, name = "reduce_logPr")
  } else {
    # minimize SS error
    sum_of_squares <- noise^-2*(expression - fit_expression)^2
    mean_squared_error <- tf$reduce_mean(sum_of_squares,
                                         axis = 0L,
                                         name = "MSE")

    loss <- tf$reduce_sum(mean_squared_error, name = "reduce_MSE")
  }
  train <- optimizer$minimize(loss = loss, name = "train")

  # train the model

  if ("noise" %in% colnames(one_timecourse)) {
    noise_entry <- matrix(
      one_timecourse$noise,
      nrow = nrow(one_timecourse),
      ncol = n_initializations
      )
  } else {
    # constant noise
    noise_entry <- matrix(
      0.2,
      nrow = nrow(one_timecourse),
      ncol = n_initializations
    )
  }

  timecourse_dict <- dict(timepts = matrix(one_timecourse$time,
                                           nrow = nrow(one_timecourse),
                                           ncol = n_initializations),
                          expression = matrix(one_timecourse$abundance,
                                              nrow = nrow(one_timecourse),
                                              ncol = n_initializations),
                          noise = noise_entry
                          )

  sess <- tf$compat$v1$Session()
  # initialize parameters
  sess$run(tf$compat$v1$global_variables_initializer())

  # keep track of initialization for error checking
  initial_vals <- lapply(parameters,
                         function(variable){
                           tibble::tibble(
                             variable = variable,
                             init_id = 1:n_initializations,
                             value = sess$run(eval(parse(text = variable)))
                           )
                         }) %>%
    dplyr::bind_rows()

  # find an NLS/MAP minima for each initialization

  past_loss <- 100000
  continue <- TRUE
  while (continue) {
    # train
    for (i in 1:100) {
      sess$run(
        train,
        feed_dict = timecourse_dict
        )
    }

    # loss for individual parameter sets
    current_losses <- if (use_prior) {
      sess$run(logPr, feed_dict = timecourse_dict)
    } else {
      sess$run(mean_squared_error, feed_dict = timecourse_dict)
    }

    if (sum(!is.nan(current_losses)) < pmin(10, n_initializations)) {
      # return NULL so the run can be re-initialized
      return(NULL)
    } else {
      valid_summed_loss <- sum(current_losses, na.rm = TRUE)

      if (verbose) {
        print(valid_summed_loss)
      }

      if (past_loss - valid_summed_loss > 0.001) {
        past_loss <- valid_summed_loss
      } else{
        continue <- FALSE
      }
    }
  }

  # summarize valid (and invalid) parameter sets

  output <- list()

  # invalid parameter set initial parameters

  if (any(is.nan(current_losses))) {
    output$invalid_timecourse_fits <- initial_vals %>%
      dplyr::filter(init_id %in% which(is.nan(current_losses))) %>%
      dplyr::mutate(tc_id = a_tc_id) %>%
      dplyr::select(tc_id, init_id, variable, value)
  } else {
    output$invalid_timecourse_fits <- data.frame()
  }

  # valid parameter set optimal parameters, fits, MSE
  valid_parameter_sets <- which(!is.nan(current_losses))

  # fit parameters

  output$parameters <- lapply(
    parameters,
    function(variable){
      tibble::tibble(variable = variable,
                     init_id = 1:n_initializations,
                     value = sess$run(eval(parse(text = variable))))
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::filter(init_id %in% valid_parameter_sets) %>%
    dplyr::mutate(tc_id = a_tc_id) %>%
    dplyr::select(tc_id, init_id, variable, value)

  output$loss <- if (use_prior) {
    tibble::tibble(tc_id = a_tc_id,
                   init_id = valid_parameter_sets,
                   loss = current_losses[valid_parameter_sets],
                   logLik = sess$run(
                     normal_logLik,
                     feed_dict = timecourse_dict)[valid_parameter_sets],
                   logPriorPr = sess$run(
                     model_log_pr,
                     feed_dict = timecourse_dict)[valid_parameter_sets])
  } else {
    tibble::tibble(tc_id = a_tc_id,
                   init_id = valid_parameter_sets,
                   loss = current_losses[valid_parameter_sets])
  }

  return(output)
}

init_sigmoid_parameters <- function(use_prior, n_initializations, prior_pars = NULL) {

  checkmate::assertLogical(use_prior, len = 1)
  checkmate::assertNumber(n_initializations, lower = 10)

  if (use_prior) {
    v_inter <- tf$Variable(
      tf$random$normal(shape(n_initializations),
                       mean = 0,
                       stddev = prior_pars["v_sd"]),
      name = "v_inter")
    t_rise <- tf$Variable(
      tf$random$gamma(shape(n_initializations),
                      alpha = prior_pars["time_shape"],
                      beta = 1 / prior_pars["time_scale"]),
      name = "t_rise")
    rate <- tf$Variable(
      tf$random$gamma(shape(n_initializations),
                      alpha = prior_pars["rate_shape"],
                      beta = 1 / prior_pars["rate_scale"]),
      name = "rate")
  } else {
    v_inter <- tf$Variable(
      tf$random$normal(shape(n_initializations),
                       mean = 0,
                       stddev = prior_pars["v_sd"]),
      name = "v_inter")
    t_rise <- tf$Variable(
      tf$random$uniform(shape(n_initializations),
                        0,
                        prior_pars["t_max"]),
      name = "t_rise")
    rate <- tf$Variable(
      tf$random$uniform(shape(n_initializations),
                        0,
                        1),
      name = "rate")
  }

  return(
    list(v_inter = v_inter,
         t_rise = t_rise,
         rate = rate)
  )
}

init_impulse_parameters <- function(use_prior, n_initializations, prior_pars = NULL, t_rise) {

  checkmate::assertLogical(use_prior, len = 1)
  checkmate::assertNumber(n_initializations, lower = 10)

  if (use_prior) {
    v_final <- tf$Variable(
      tf$random$normal(shape(n_initializations),
                       mean = 0,
                       stddev = prior_pars["v_sd"]),
      name = "v_final")
    t_diff <- tf$Variable(
      tf$random$gamma(shape(n_initializations),
                      alpha = prior_pars["time_shape"],
                      beta = 1 / prior_pars["time_scale"]),
      name = "t_diff")
  } else {
    v_final <- tf$Variable(
      tf$random$normal(shape(n_initializations),
                       mean = 0,
                       stddev = prior_pars["v_sd"]),
      name = "v_final")
    t_diff <- tf$Variable(
      tf$random$uniform(shape(n_initializations),
                        0,
                        prior_pars["t_max"]),
      name = "t_diff")
  }
  t_fall <- tf$add(t_rise, t_diff, name = "t_final")

  return(
    list(
      v_final = v_final,
      t_diff = t_diff,
      t_fall = t_fall
    )
  )
}

#' Reduce to best timecourse parameters
#'
#' Across multiple fits of a timecourse summarize the best fitting timecourse in terms of least-squares error as well as by lowest absolute V within a tolerance of the least-squares set.
#'
#' @param timecourse_list List output from \code{\link{estimate_timecourse_params_tf}}
#' @param reduction_type How to choose the best parameter set, options are:
#' \itemize{
#'  \item{\code{loss-min}: lowest loss function},
#'  \item{\code{loss-small-v-small}: loss within \code{sufficiency_tolerance} of minimum loss and then minimize absolute sum of \eqn{v_{inter}} and \eqn{v_{final}} (useful primarily when not using priors).}
#' }
#' @param sufficiency_tolerance All timecourses within 1 + sufficiency_tolerance best fitting parameter set are deemed sufficient
#'
#' @return a list containing top parameter set and losses
#'
#' @export
reduce_best_timecourse_params <- function(timecourse_list,
                                          reduction_type = "loss-min",
                                          sufficiency_tolerance = 0.05) {

  stopifnot(class(timecourse_list) == "list")
  stopifnot(names(timecourse_list) == c("invalid_timecourse_fits",
                                        "parameters", "loss"))
  if ("model" %in% colnames(timecourse_list$loss) ||
      "model" %in% timecourse_list$parameters) {
    stop ("\"model\" cannot be included in tables of timecourse list")
  }

  stopifnot(class(sufficiency_tolerance) %in% c("numeric", "integer"),
            length(sufficiency_tolerance) == 1,
            sufficiency_tolerance >= 0)

  stopifnot(class(reduction_type) == "character",
            length(reduction_type) == 1)
  valid_reduction_types <- c("loss-min", "loss-small-v-small")
  if (!(reduction_type %in% valid_reduction_types)) {
    stop (reduction_type, " is an invalid \"reduction_type\",
          valid types are: ", paste(valid_reduction_types, collapse = ", "))
  }

  good_inits <- timecourse_list$loss %>%
    dplyr::group_by(tc_id) %>%
    dplyr::arrange(loss) %>%
    dplyr::mutate(n_near_min = sum(loss - min(loss) < sufficiency_tolerance),
                  all_valid = dplyr::n())

  good_init_parameters <- good_inits %>%
    # select all parameter sets within 0.05 of "best" fitting parameter set
    dplyr::filter(loss - min(loss) < sufficiency_tolerance) %>%
    # add all parameter sets within ss tolerance
    dplyr::left_join(timecourse_list$parameters, by = c("tc_id", "init_id"))

  # range of parameter values for "good parameter sets"
  parameter_range <- good_init_parameters %>%
    dplyr::group_by(tc_id, variable) %>%
    dplyr::summarize(
      min_value = min(value),
      max_value = max(value),
      .groups = "drop"
      ) %>%
    dplyr::ungroup()

  if (reduction_type == "loss-min") {

    # lowest loss
    best_inits <- good_inits %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

  } else if (reduction_type == "loss-small-v-small") {

    # min absolute sum of v within a tolerance of best fit
    parameter_set_v_abs_sum <- good_init_parameters %>%
      dplyr::group_by(tc_id, init_id) %>%
      dplyr::summarize(
        v_abs_sum = sum(abs(value[variable %in% c("v_inter", "v_final")])),
        .groups = "drop"
        )

    parameter_set_v_abs_sum <- parameter_set_v_abs_sum %>%
      dplyr::group_by(tc_id) %>%
      dplyr::arrange(v_abs_sum) %>%
      dplyr::slice(1)

    best_inits <- good_inits %>%
      dplyr::inner_join(parameter_set_v_abs_sum, by = c("tc_id", "init_id")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(top_hit_type = "low loss, low absolute v")

  } else {
    stop ("\"", reduction_type, "\" logic not defined")
  }

  top_timecourse_fits <- list()
  top_timecourse_fits$parameters <- timecourse_list$parameters %>%
    dplyr::semi_join(best_inits, by = c("tc_id", "init_id")) %>%
    dplyr::left_join(parameter_range, by = c("tc_id", "variable"))
  top_timecourse_fits$loss <- best_inits

  top_timecourse_fits
}

#' Fit timecourse
#'
#' Fit the parameters of a sigmoid or impulse model to a set of time points.
#'
#' @param timecourse_parameters a one row tibble containing each kinetic parameter as a separate column.
#' @param timepts a numeric vector of timepoints to fit
#' @param model sigmoid or impulse
#' @param fit.label fitted values variable name
#'
#' @examples
#' timecourse_parameters <- tibble::tibble(t_rise = 25, rate = 0.25, v_inter = 3,
#'                                        v_final = -3, t_fall = 45)
#' timecourse_parameters <- tibble::tibble(t_rise = 45, rate = 1, v_inter = 3)
#' fit_timecourse(timecourse_parameters, model = "sigmoid")
#'
#' @export
fit_timecourse <- function (
    timecourse_parameters,
    timepts = seq(0, 90, by = 1),
    model = "sigmoid",
    fit.label = "fit"
    ) {

  stopifnot("data.frame" %in% class(timecourse_parameters),
            nrow(timecourse_parameters) == 1)
  stopifnot(all(class(timepts) %in% c("numeric", "integer")),
            length(timepts) > 0)

  stopifnot(class(model) == "character", length(model) == 1)
  if (model == "sigmoid") {
    stopifnot(all(c("v_inter", "t_rise", "rate") %in%
                    colnames(timecourse_parameters)))
  } else if (model == "impulse") {
    stopifnot(all(c("v_inter", "t_rise", "rate", "v_final", "t_fall") %in%
                    colnames(timecourse_parameters)))
  } else {
    stop(model, ' is not a valid option for "model",
         use "sigmoid" or "impulse"')
  }
  stopifnot(class(fit.label) == "character", length(fit.label) == 1)

  # combine parameters + times

  eval_times <- timecourse_parameters %>%
    dplyr::mutate(time = purrr::map(1, function(x){timepts})) %>%
    tidyr::unnest(time)

  eqtn <- switch(model,
                 "sigmoid" = "v_inter*(1/(1 + exp(-1*rate*(time - t_rise))))",
                 "impulse" = "(1/(1 + exp(-1*rate*(time - t_rise)))) *
                 (v_final + (v_inter - v_final)*
                   (1/(1 + exp(rate*(time - t_fall)))))")

  eval_times$fit <- eval(parse(text = eqtn), eval_times)

  fitted_timecourse <- eval_times %>%
    dplyr::select(time, !!fit.label := fit)

  if ("tzero_offset" %in% colnames(timecourse_parameters)) {
    # adjust intercept if an offset is provided
    fitted_timecourse <- fitted_timecourse %>%
      dplyr::mutate(!!fit.label := !!rlang::sym(fit.label) + timecourse_parameters$tzero_offset[[1]])
  }

  return(fitted_timecourse)
}
