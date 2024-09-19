#' @title Combine the Signal from 5hmC and 5mC
#'
#' @description Merges methylation signal for both "h" and "m". This is useful
#'  when comparing ONT data to WGBS data, as WGBS doesn't discern between 5mC
#'  and 5hmC (unlike modifided base called ONT data).
#'
#' @inheritParams comparative_bedmethyl_format
#'
#' @return A comparitive bedmethyl file with merged 5mC and 5hmC signal
#'
#' @examples
#' # Read in a comparative bedmethyl file
#' bedmethyl <- read_comparison_bedmethyl("path/to/bedmethyl.bed")
#'
#' head(bedmethyl, 4)
#'       chr start   end mark_name ONT_read_depth ONT_percent_methylation
#'    <char> <num> <num>    <char>          <num>                   <num>
#' 1:   chr1     0   200         h              6                       5
#' 2:   chr1     0   200         m              6                       5
#' 3:   chr1   200   400         h             99                      10
#' 4:   chr1   200   400         m             99                      20
#'    BS_read_depth BS_percent_methylation
#'            <num>                  <num>
#' 1:             0                      0
#' 2:             0                      0
#' 3:             0                      0
#' 4:             0                      0
#'
#' # Combine signal
#' bedmethyl <- combine_5mc_5hmc(bedmethyl)
#'
#' head(bedmethyl, 2)
#'       chr start   end mark_name ONT_read_depth ONT_percent_methylation
#'    <char> <num> <num>    <char>          <num>                   <num>
#' 1:   chr1     0   200         m              6                      10
#' 2:   chr1   200   400         m             99                      30
#'    BS_read_depth BS_percent_methylation
#'            <num>                  <num>
#' 1:             0                      0
#' 2:             0                      0
#' @export
combine_5mc_5hmc <- function(comparison_bedmethyl) {
  methylation_5mc <- dplyr::filter(
    comparison_bedmethyl,
    mark_name == "m"
  )
  methylation_5hmc <- dplyr::filter(
    comparison_bedmethyl,
    mark_name == "h"
  ) |>
    dplyr::mutate(
      mark_name = "m"
    )

  # ONT data (as called by `modkit pileup`) displays the read depth for each
  # site and the percentage of reads that are 5hmC and 5mC. To combine them
  # we therefore just need to add these percentages (as they can never add to
  # a value above 100%)
  comparison_bedmethyl <- dplyr::full_join(
    methylation_5mc,
    methylation_5hmc,
    by = c(
      "chr",
      "start",
      "end",
      "mark_name",
      "ONT_read_depth",
      "BS_read_depth",
      "BS_percent_methylation"
    )
  ) |>
    dplyr::mutate(
      ONT_percent_methylation = ONT_percent_methylation.x +
        ONT_percent_methylation.y
    ) |>
    dplyr::select(
      c(
        "chr",
        "start",
        "end",
        "mark_name",
        "ONT_read_depth",
        "ONT_percent_methylation",
        "BS_read_depth",
        "BS_percent_methylation"
      )
    )

  return(comparison_bedmethyl)
}
