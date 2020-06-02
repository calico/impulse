#' Kinetics plotting
#'
#' @param augmented_timecourses Tibble with one row per timecourse and containing measurements &/or fitted kinetics to plot.
#' @param saturation Degree of saturation to use for shading sigmoidal curves.
#' @param max_time Maximum time to display
#' @param fit_timepoints Number of timepoints to use for generating fitted timecourses.
#'
#' @return a ggplot
#'
#' @export
kinetics_plotting <- function(augmented_timecourses,
                              saturation = 0.9,
                              max_time,
                              fit_timepoints = 100) {

  stopifnot("tbl_df" %in% class(augmented_timecourses))
  stopifnot("tc_id" %in% colnames(augmented_timecourses))

  reserved_variables <- c("measurements", "fitted_kinetics")
  present_variables <- intersect(reserved_variables,
                                 colnames(augmented_timecourses))
  if (length(present_variables) == 0) {
    stop ("no aesthetic variables present; aesthetic variables are: ",
          paste(reserved_variables, collapse = ", "))
  }

  stopifnot(class(saturation) == "numeric",
            length(saturation) == 1,
            saturation > 0.5,
            saturation < 1)
  stopifnot(class(max_time) %in% c("numeric", "integer"),
            length(max_time) == 1,
            max_time > 0)
  stopifnot(class(fit_timepoints) %in% c("numeric", "integer"),
            length(fit_timepoints) == 1,
            fit_timepoints > 0)

  if ("measurements" %in% present_variables) {
    measurements <- augmented_timecourses %>% tidyr::unnest(measurements)
    stopifnot(all(c("tc_id", "time", "abundance") %in% colnames(measurements)))

    kinetics_plot <- ggplot(measurements) +
      geom_point(aes(x = time, y = abundance)) +
      facet_wrap(~ tc_id)

  } else {
    # no measurements provide
    kinetics_plot <- ggplot(tibble::tibble(tc_id = NA_integer_[-1],
                                           time = 0[-1],
                                           abundance = 0[-1])) +
      geom_point(aes(x = time, y = abundance)) +
      facet_wrap(~ tc_id)
  }

  if ("fitted_kinetics" %in% present_variables) {

    fitted_kinetics <- augmented_timecourses %>%
      dplyr::select(tc_id, best_model, fitted_kinetics) %>%
      tidyr::unnest_legacy(fitted_kinetics)

    timepoints <- seq(from = 0, max_time, length.out = fit_timepoints)
    fitted_values <- fitted_kinetics %>%
      tidyr::nest_legacy(rate, t_rise, t_fall, v_inter, v_final,
                  .key = "parameters") %>%
      dplyr::mutate(fitted_timecourses = purrr::map2(parameters, model,
                                                     fit_timecourse,
                                                     timepts = timepoints)) %>%
      tidyr::unnest_legacy(fitted_timecourses) %>%
      dplyr::mutate(best_model = ifelse(model == best_model, TRUE, FALSE))

    # calculate saturation interval

    kinetic_intervals <- kinetic_aesthetics(
      fitted_kinetics %>%
        dplyr::filter(model == best_model),
      fitted_values,
      saturation = saturation
      )

    kinetics_plot <- kinetics_plot +
      geom_rect(data = kinetic_intervals$asympote_aesthetics_df,
                aes(xmin = t_saturation_start,
                    xmax = t_saturation_end),
                ymin = -Inf, ymax = Inf, alpha = 0.25) +
      geom_hline(yintercept = 0, color = "BLACK", size = 1) +
      geom_path(data = fitted_values, aes(x = time, y = fit,
                                          color = model, group = model,
                                          linetype = best_model)) +
      # timing
      geom_segment(data = kinetic_intervals$time_aesthetics_df,
                   aes(x = time, y = fit - 0.25,
                       xend = time, yend = fit + 0.25),
                   color = "limegreen", size = 2) +
      geom_text(data = kinetic_intervals$time_aesthetics_df %>%
                  dplyr::filter(par == "t_rise"), aes(x = time + 2, y = fit),
                label = "t[rise]", parse = TRUE, hjust = 0, size = 4) +
      geom_text(data = kinetic_intervals$time_aesthetics_df %>%
                  dplyr::filter(par == "t_fall"), aes(x = time + 2, y = fit),
                label = "t[fall]", parse = TRUE, hjust = 0, size = 4) +
      # assymptotes
      geom_segment(data = kinetic_intervals$asympote_aesthetics_df,
                   aes(x = t_saturation_end - 5, y = assymp,
                       xend = t_saturation_end + 5, yend = assymp),
                   color = "limegreen", size = 2) +
      geom_text(data = kinetic_intervals$asympote_aesthetics_df %>%
                  dplyr::filter(assymp_type == "v_inter"),
                aes(x = t_saturation_end,
                    y = assymp + ifelse(assymp >= 0, 0.1, -0.1)),
                label = "v[inter]", parse = TRUE, size = 4) +
      geom_text(data = kinetic_intervals$asympote_aesthetics_df %>%
                  dplyr::filter(assymp_type == "v_final"),
                aes(x = t_saturation_end,
                    y = assymp + ifelse(assymp >= 0, 0.1, -0.1)),
                label = "v[final]", parse = TRUE, size = 4) +
      scale_linetype_manual(values = c("TRUE" = 1, "FALSE" = 2)) +
      coord_cartesian(xlim = c(0, max_time))

  }

  kinetics_plot
}

kinetic_aesthetics <- function(fitted_kinetics,
                               fitted_values,
                               saturation = 0.9) {

  # find fitted kinetic value at timing coefficient
  time_aesthetics_df <- fitted_kinetics %>%
    dplyr::select(tc_id, model, t_rise, t_fall) %>%
    tidyr::gather(par, value, -tc_id, -model) %>%
    dplyr::filter(!is.na(value)) %>%
    dplyr::left_join(fitted_values %>%
                       dplyr::select(time, fit, tc_id, model),
                     by = c("tc_id", "model")) %>%
    dplyr::group_by(tc_id, par) %>%
    dplyr::arrange(abs(time - value)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  # find satuation time for assymptote
  asympote_aesthetics_df <- fitted_kinetics %>%
    dplyr::select(tc_id, model, rate, t_rise, t_fall, v_inter, v_final) %>%
    tidyr::gather(par, value, -tc_id, -model, -rate, -v_inter, -v_final) %>%
    dplyr::filter(!is.na(value)) %>%
    # match t_rise to v_inter and t_fall to v_final
    dplyr::mutate(assymp = dplyr::case_when(par == "t_rise" ~ v_inter,
                                            par == "t_fall" ~ v_final)) %>%
    dplyr::mutate(assymp_type = dplyr::case_when(
      par == "t_rise" ~ "v_inter",
      par == "t_fall" ~ "v_final")) %>%
    # derive saturation time for x-axis
    dplyr::mutate(t_saturation_start = saturation_time(1 - saturation,
                                                       value,
                                                       rate),
                  t_saturation_end = saturation_time(saturation,
                                                     value,
                                                     rate)) %>%
    dplyr::select(tc_id, model, assymp_type, assymp,
                  t_saturation_start, t_saturation_end)

  list(time_aesthetics_df = time_aesthetics_df,
       asympote_aesthetics_df = asympote_aesthetics_df)
}

saturation_time <- function(saturation, c_time, rate) {

  c_time + log(saturation / (1 - saturation)) / rate

}
