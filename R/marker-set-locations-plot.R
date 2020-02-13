#' Plot genomic locations of different marker sets for comparison
#'
#' Pass this thing a tibble of markers and it will make a plot with
#' the marker names in the order they are given, with lines to
#' the genomic coordinate relative to the reference used. Assumes that all
#' markers are on the same chromosome.
#'
#' Function returns a ggplot object.  To change colors of the marker sets,
#' add a new discrete fill scale to it, as shown in the examples. Be sure to
#' add \code{na.value = "transparent"}  and set breaks to be the reverse of the
#' names of columns 2 and up (see examples).
#'
#' @param M the tibble of markers.  Must be in a particular format: 1) First column
#' is numeric and gives the genomic coordinate of the variant/marker in each row.
#' The reference genome used determines the x-coordinate of the plot. The name of this
#' column must be "position". 2) the second column must be named pos_predicted and be a logical vector that
#' says whether the position of the marker in the first column is predicted or not.
#' If not, set them all to FALSE.
#' 3) Each successive
#' column is used to store the markers names
#' in a particular marker set.  The title of the column should be the name of the marker
#' set.  The name of each marker in the marker set should appear in the row which is at the
#' appropriate location at the corresponding position in the coordinate column.  4) There must
#' be no other columns in the tibble.
#' 5) The second column could be a character vector
#' that is the chrom-name:position for the reference used.  6) The order of the marker sets/coordinates
#' in columns 2 and up
#' determines the order in which they will be placed going from the bottom to the top of
#' the plot. 7) Note that the position of some markers may not be known on the
#' reference coordinates, typically because they are from a different assembly and they
#' do not align to the reference assembly.  Nonetheless, their order within the
#' rows of markers should be reasonably well defined.  Put them in in the correct order
#' and leave their position row as NA.  Note that the "outermost" markers should have
#' non-NA positions.
#' @param xlabel  The desired label beneath the x-axis.  Typically will refer to what the
#' coordinate reference is, such as, "Position on chr28 from Otsh_v1.0". If left empty,
#' it will default to "position"
#' @param linehang a 2-vector (or recycled if shorter) that gives the overhang amount of the
#' position line to the left of the leftmost markers (element 1) and to the right of the rightmost
#' marker, expressed as a proportion of the total distance between leftmost and rightmost marker.
#' @param boxhang a 2-vector giving the fraction of the distance between outermost markers that
#' we want the box edges to be at on the left and right.  Negative moves them in toward the center
#' and positive moves them out.
#' @param m_name_upnudge absolute amount by which to nudge the marker names upward within their boxes.
#' @param m_name_textsize size of text for the marker names.  This should be fiddled with *after*
#' the output size has been set, to make the marker names fit in their boxes.
#' @param kink_rel the height of the kinks (and of the ticks on the position line) as a fraction
#' of the total height of all the marker boxes.
#' @param conn_line_size the size (thickness) of the lines/kinks/ticks connecting marker boxes to their positions.
#' @export
#' @examples
#' # mykiss markers plot
#' # none of these positions are predicted, but this is for testing.
#' M <- readr::read_csv("inst/extdata/mykiss-greb1l-snp-summary.csv") %>%
#'        mutate(pos_predicted = rep(c(FALSE, FALSE, TRUE, FALSE), length.out = dplyr::n())) %>%
#'        select(position, pos_predicted, everything())
#'
#' g1 <- marker_set_locations_plot(M)
#' g2 <- marker_set_locations_plot(M) +
#'   scale_fill_brewer(
#'     palette = "Set2",
#'     na.value = "transparent",
#'     breaks = rev(names(M)[-(1:2)])
#'    )
#'
#' # for a final figure, choose your output size and format
#' # and work from that to set the m_name_textsize
#' \notrun{
#' g3 <- marker_set_locations_plot(M,
#'         m_name_upnudge = 0.007,
#'         m_name_textsize = 3.0
#' ) +
#'   scale_fill_brewer(
#'     palette = "Set2",
#'     na.value = "transparent",
#'     breaks = rev(names(M)[-1])
#'    )
#'  ggsave(g3,
#'  filename = "msl-plot.pdf",
#'  width = 8, height = 10
#'  )
#' }
marker_set_locations_plot <- function(M,
                                      xlabel = "position",
                                      linehang = 0.07,
                                      boxhang = -0.1,
                                      m_name_upnudge = 0.01,
                                      m_name_textsize = 3.0,
                                      kink_rel = 0.08,
                                      conn_line_size = 0.3
                                      ) {

  NMB <- nrow(M)  # number of marker boxes

  # figure out x positions for things.
  # the line at the bottom to place positions on
  if (length(linehang) < 2) {
    linehang <- rep(linehang, times = 2)
  }
  mextent <- range(M$position, na.rm = TRUE)
  mlo <- mextent[1]
  mhi <- mextent[2]
  mdiff <- mhi - mlo
  linelo <- mlo - linehang[1] * mdiff
  linehi <- mhi + linehang[2] * mdiff

  # the left points and widths of the marker boxes
  if (length(boxhang) < 2) {
    boxhang <- rep(boxhang, times = 2)
  }
  boxlo <- mlo - boxhang[1] * mdiff
  boxhi <- mhi + boxhang[2] * mdiff
  boxwidth <- (boxhi - boxlo) / NMB

  # Now, prep the tibble...
  # get levels for the marker sets:
  ms_levs <- names(M)[-1]

  # put the left and right end points of the boxes on
  M2 <- M %>%
    dplyr::mutate(box_lefts = boxlo + (0:(dplyr::n() - 1)) * boxwidth,
                  box_rights = box_lefts + boxwidth,
                  box_mids = box_lefts + 0.5 * boxwidth
                  ) %>%
    tidyr::gather(key = "marker_set", value = "marker_name",
                  -position, -pos_predicted, -box_lefts, -box_rights, -box_mids) %>%
    dplyr::mutate(marker_set = factor(marker_set, levels = ms_levs)) %>%
    dplyr::group_by(marker_set) %>%
    dplyr::mutate(max_name_strwidth = max(strwidth(marker_name, units = "figure"), na.rm = TRUE))


  # Now, set the heights of the boxes according to the max strwidths + some padding
  M3 <- M2 %>%
    dplyr::summarise(height = max_name_strwidth[1]) %>%
    dplyr::mutate(
                  box_bottoms = 1.5 * max(height) + cumsum(dplyr::lag(height, default = 0)),
                  box_tops = box_bottoms + height
                  ) %>%
    dplyr::select(marker_set, box_bottoms, box_tops) %>%
    dplyr::left_join(M2, ., by = "marker_set") %>%
    dplyr::mutate(`Marker Set` = marker_set) %>%
    dplyr::ungroup()

  M3$`Marker Set`[is.na(M3$marker_name)] <- NA

  # get the total height of all boxes
  total_boxes_height <- max(M3$box_tops) - min(M3$box_bottoms)
  kink_ht <- total_boxes_height * kink_rel

  M4 <- M3 %>%
    dplyr::mutate(tick_top = kink_ht * 0.5,
           tick_bottom = -kink_ht * 0.5,
           kink_top = min(box_bottoms),
           kink_bottom = kink_top - kink_ht
           )


  ggplot(M4) +
    geom_rect(
              aes(xmin = box_lefts,
                  xmax = box_rights,
                  ymin = box_bottoms,
                  ymax = box_tops,
                  fill = `Marker Set`),
              colour = "black"
              ) +
    geom_text(
              aes(x = box_mids, y = box_bottoms, label = marker_name),
              hjust = 0,
              nudge_y = m_name_upnudge,
              angle = 90,
              size = m_name_textsize
              ) +
    annotate("segment", x = linelo, xend = linehi, y = 0, yend = 0) + # position line
    geom_segment(aes(x = position, xend = position, y = tick_bottom, yend = tick_top),
                 size = conn_line_size) +  # tick on the position line
    geom_segment(aes(x = box_mids, xend = box_mids, y = kink_bottom, yend = kink_top, colour = pos_predicted),
                 size = conn_line_size) +  # kinks below bottom marker set boxes
    geom_segment(aes(x = position, xend = box_mids, y = tick_top, yend = kink_bottom, colour = pos_predicted),
                 size = conn_line_size) +  # connecting lines
    scale_fill_discrete(na.value = "transparent",
                        breaks = rev(ms_levs)) +
    theme_bw() +
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          panel.grid = element_blank(),
          axis.ticks.y = element_blank()) +
    xlab(xlabel) +
    scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "black"))
}
