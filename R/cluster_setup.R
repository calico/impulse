#' Setup timecourse fit jobs
#'
#' @examples
#' timecourse_fit_job_submission(expression_data, njobs = 50, remote_path = "sean@sci-vm-005:/home/sean/data/McIsaac_timecourses/timecourse_data.Rds")
#'
#' @export
timecourse_fit_job_submission <- function(expression_data, njobs = 200, remote_path) {

  stopifnot(all(colnames(expression_data) %in% c("TF", "strain", "date", "restriction", "mechanism", "time", "GeneName", "log2_fc")))

  timecourses <- expression_data %>%
    dplyr::select(TF, strain, date, restriction, mechanism, time, GeneName, log2_fc) %>%
    dplyr::group_by(TF, strain, date, restriction, mechanism, GeneName) %>%
    dplyr::filter(!(all(log2_fc == 0)))

  timecourse_labels <- timecourses %>%
    dplyr::ungroup() %>%
    dplyr::distinct(TF, strain, date, restriction, mechanism, GeneName) %>%
    dplyr::mutate(tc_id = 1:n()) %>%
    dplyr::mutate(job_id = c(rep(1:njobs, times = floor(nrow(.)/njobs)), 1:(nrow(.) %% njobs)))

  timecourses <- timecourses %>%
    dplyr::left_join(timecourse_labels, by = c("TF", "strain", "date", "restriction", "mechanism", "GeneName")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(time = as.numeric(as.character(time)))

  tmp_file = file.path("/tmp", "timecourse_data.Rds")
  saveRDS(timecourses, file = tmp_file)
  cat(paste0('rsync -avP ', tmp_file, ' ', remote_path))
}
