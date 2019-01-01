

#' expand SAM flag into 12 columns
#'
#' @param D the tibble of SAM rows
#' @param column the name of the column holding the SAM flag
#' @export
expand_sam_flag <- function(D, column) {
  flags <- D[[column]]
  m <- simplify2array(lapply(flags,function(x){ as.logical(intToBits(x))}))[1:12,] %>%
      t(.) %>%
      tibble::as_tibble() %>%
      setNames(c("flag_paired",
                 "flag_paired_proper",
                 "flag_unmapped",
                 "flag_mate_unmapped",
                 "flag_read_reverse_strand",
                 "flag_mate_reverse_strand",
                 "flag_first_in_pair",
                 "flag_second_in_pair",
                 "flag_not_primary_align",
                 "flag_fails_platform_check",
                 "flag_pcr_duplicate",
                 "flag_is_supplemental"
      ))

  dplyr::bind_cols(D, m)

}
