#' Impulse Sigmoid Comparison
#'
#' @param timecourse_losses table containing tc_id, model and logLik
#' @param fdr_cutoff cutoff for accepting an impulse over a sigmoid fit
#'
#' @export
impulse_sigmoid_comparison <- function(timecourse_losses, fdr_cutoff = 0.001) {

  stopifnot("data.frame" %in% class(timecourse_losses))
  required_colnames <- c("tc_id", "model", "logLik")
  missing_required_colnames <- setdiff(required_colnames, colnames(timecourse_losses))
  if (length(missing_required_colnames) != 0) {
    stop ("\"timecourse_losses\" is missing required variables: ", paste(missing_required_colnames, collapse = ", "))
  }

  stopifnot(class(fdr_cutoff) == "numeric", fdr_cutoff > 0, fdr_cutoff < 1)

  spread_models <- timecourse_losses %>%
    dplyr::select(tc_id, model, logLik) %>%
    tidyr::spread(model, logLik)

  required_vars_spread_models <- c("tc_id", "sigmoid", "impulse")
  missing_required_vars_spread_models <- setdiff(required_vars_spread_models, colnames(spread_models))
  if (length(missing_required_vars_spread_models) != 0) {
    stop ("\"timecourse_losses\" is missing required models: ", paste(missing_required_vars_spread_models, collapse = ", "))
  }

  incomplete_comparisons <- spread_models %>%
    dplyr::filter(is.na(sigmoid) | is.na(impulse))

  if (nrow(incomplete_comparisons) != 0) {
    warning (nrow(incomplete_comparisons), " tc_ids do not contain both a sigmoid and impulse model")
  }

  spread_models %>%
    dplyr::mutate(logLik_diff = impulse - sigmoid,
                  logLik_diff = pmax(logLik_diff, 0),
                  # apply likelihood ratio test
                  model_pchisq = stats::pchisq(logLik_diff, df = 2, lower.tail = FALSE),
                  # correct using Benjamini-Hochberg (using this instead of Storey because there is an enrichment of p-values near 1 due to
                  # timecourses with strong prior constraints which were fit much better with a sigmoid than an impulse
                  model_qchisq = stats::p.adjust(model_pchisq, method = "BH"),
                  best_model = dplyr::case_when(model_qchisq < 0.001 ~ "impulse",
                                                TRUE ~ "sigmoid")) %>%
    dplyr::select(tc_id, logLik_diff, model_pchisq, model_qchisq, best_model)
}
