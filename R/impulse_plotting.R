sigmoid_impulse_plot <- function(timecourse_parameters) {

  fit_timecourse(timecourse_parameters, model = "sigmoid") %>%
    dplyr::rename(sigmoid = fit) %>%
    dplyr::left_join(fit_timecourse(timecourse_parameters, model = "impulse") %>%
                       dplyr::rename(impulse = fit), by = "time") %>%
    tidyr::gather(eqtn, level, -time) %>%
    ggplot(aes(x = time, y = level, color = eqtn)) + geom_path(size = 2) +
    theme_bw()

}

#' Gene x TF Timecourse
#'
#' @inheritParams gene_tf_timecourse
#' @param spread_kinetics sigmoids and impulse kinetic summaries for each timecourse
#' @param quant_var a variable to use for fold changes
#'
#' @return a ggplot2 plot
#'
#' @export
gene_tf_fits <- function(expression_tall_data,
                         spread_kinetics,
                         query_TF,
                         query_gene,
                         quant_var = "shrunken",
                         time_trans = "identity") {

  # summary of all measurements

  expression_measurements <- expression_tall_data %>%
    dplyr::filter(TF == query_TF, gene == query_gene) %>%
    dplyr::mutate(tf_label = paste(strain, date, restriction, mechanism, sep = "-"),
                  time = as.numeric(as.character(time))) %>%
    dplyr::rename(log2fc = !!rlang::sym(quant_var)) %>%
    dplyr::select(tc_id, tf_label, time, log2fc)

  # find matching parametric fits

  reduced_kinetics <- spread_kinetics %>%
    dplyr::inner_join(expression_measurements %>%
                        dplyr::distinct(tc_id, tf_label), by = "tc_id") %>%
    dplyr::mutate(kinetics = ifelse(is_impulse, "impulse", "sigmoid"))

  if (nrow(reduced_kinetics) != 0) {

    timecourse_fits <- reduced_kinetics %>%
      plyr::dlply(.variables = c("tc_id", "model")) %>%
      lapply(function(x){
        dynamicyeast::fit_timecourse(x, timepts = seq(0, max(expression_measurements$time)), model = x$model) %>%
          dplyr::mutate(tc_id = x$tc_id, tf_label = x$tf_label, model = x$model, is_model = x$model == x$kinetics)
      }) %>%
      dplyr::bind_rows()

    reduced_kinetics_is_model <- reduced_kinetics %>%
      dplyr::filter(is_impulse & model == "impulse" | !is_impulse & model == "sigmoid")

    time_aesthetics_df <- reduced_kinetics_is_model %>%
      dplyr::select(tc_id, model, tf_label, t_rise, t_fall) %>%
      tidyr::gather(par, value, -tc_id, -tf_label, -model) %>%
      dplyr::filter(!is.na(value)) %>%
      dplyr::left_join(timecourse_fits %>%
                         dplyr::select(time, fit, tc_id, model), by = c("tc_id", "model")) %>%
      dplyr::group_by(tc_id, par) %>%
      dplyr::arrange(abs(time - value)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

    asympote_aesthetics_df <- reduced_kinetics_is_model %>%
      dplyr::select(tc_id, tf_label, rate, t_rise, t_fall) %>%
      tidyr::gather(par, value, -tc_id, -tf_label, -rate) %>%
      dplyr::filter(!is.na(value)) %>%
      # match t_rise to v_inter and t_fall to v_final
      dplyr::left_join(reduced_kinetics_is_model %>%
                         dplyr::select(tc_id, v_inter, v_final),
                       by = "tc_id") %>%
      dplyr::mutate(assymp = dplyr::case_when(par == "t_rise" ~ v_inter,
                                              par == "t_fall" ~ v_final)) %>%
      dplyr::mutate(assymp_type = dplyr::case_when(par == "t_rise" ~ "v_inter",
                                                   par == "t_fall" ~ "v_final")) %>%
      # derive saturation time for x-axis
      dplyr::mutate(t_saturation = saturation_time(0.9, value, rate)) %>%
      dplyr::select(tc_id, tf_label, assymp_type, t_saturation, assymp)
  }

  # plot experimental data

  measured_data_plot <- ggplot(expression_measurements, aes(x = time))

  if (nrow(reduced_kinetics) != 0) {
    measured_data_plot <- measured_data_plot +
      geom_hline(yintercept = 0, color = "BLACK", size = 1) +
      geom_path(data = timecourse_fits, aes(y = fit, color = model, group = paste0(tc_id, model), linetype = is_model), size = 1) +
      # timing
      geom_segment(data = time_aesthetics_df, aes(x = time, y = fit - 0.25, xend = time, yend = fit + 0.25), color = "limegreen", size = 2) +
      geom_text(data = time_aesthetics_df %>% dplyr::filter(par == "t_rise"), aes(x = time + 2, y = fit), label = "t[rise]", parse = TRUE, hjust = 0, size = 7) +
      geom_text(data = time_aesthetics_df %>% dplyr::filter(par == "t_fall"), aes(x = time + 2, y = fit), label = "t[fall]", parse = TRUE, hjust = 0, size = 7) +
      # assymptotes
      geom_segment(data = asympote_aesthetics_df, aes(x = t_saturation - 5, y = assymp, xend = t_saturation + 5, yend = assymp), color = "limegreen", size = 2) +
      geom_text(data = asympote_aesthetics_df %>% dplyr::filter(assymp_type == "v_inter"), aes(x = t_saturation, y = assymp + ifelse(assymp >= 0, 0.1, -0.1)), label = "v[inter]", parse = TRUE, size = 7) +
      geom_text(data = asympote_aesthetics_df %>% dplyr::filter(assymp_type == "v_final"), aes(x = t_saturation, y = assymp + ifelse(assymp >= 0, 0.1, -0.1)), label = "v[final]", parse = TRUE, size = 7)
  }

  measured_data_plot +
    # raw points
    geom_point(aes(y = log2fc), size = 2) +
    # formatting
    facet_grid(tf_label ~ .) +
    scale_color_manual("Model", values = c("impulse" = "tomato", "sigmoid" = "dodgerblue")) +
    scale_linetype_manual("Best Model", values = c("TRUE" = 1, "FALSE" = 3)) +
    expand_limits(y = c(-2,2)) +
    scale_y_continuous(expression(log[2] ~ "expression")) +
    scale_x_continuous("Time (minutes)", trans = time_trans, breaks = unique(expression_measurements$time)) +
    theme_bw() +
    theme(text = element_text(size = 20), axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom")
}

saturation_time <- function(saturation, c_time, rate) {

  c_time + log(saturation / (1 - saturation)) / rate

}
