



#' make a raster plot of allelic types on haplotypes
#'
#' This is pretty simple, but it is a standard operation,
#' and a hassle to hand code over and over again (especially
#' getting the lines that show the positions of the markers).
#' @param D a tibble with a column \code{POS} that is the positions of the
#' markers, a column \code{hnames} which are the names of the haplotypes, and
#' a column \code{atypes} which are the types to be plotted as different colors.
#' It might also have a column \code{maxGP} which holds the genotype probability
#' of the called allele.
#' @param h_ord if not NULL (the default) this is a vecdtor giving the
#' order (from top to bottom of
#' the plot) that you want these to come out.  This will also filter the data set
#' to include only the hnames in h_ord.  If this is NULL, then the order
#' of haplotypes in the data set will be used.
#' @param pos_annot a tibble with two columns: pname and pos that give the name of the
#' label you want on the genome position bar, and the position you want it at.
#' @param fcolors a named vector to be passed as the values argument to scale_fill_manual. Note that
#' if you are including more categories in annotation_columns, then you have to include
#' them, too, for example: fcolors <- c(S = "yellow", F = "blue", spring = "yellow", fall = "blue", Sacto = "red", Klamath = "green")
#' @param highlight_pos a tibble with a column named "POS" of SNP positions and column named "name"
#' which holds the name for SNPs that should be highlighted by surrounding
#' the columns by a dotted line.
#' @param pos_bar_height what fraction of the vertical plot area should we find the position
#' line above the plot area?
#' @param kink_frac the fraction of the pos_bar_height at which the marker lines should kink.
#' @param tick_frac the fraction of pos_bar_height that will be half the tick mark height.
#' @param expand_top the fraction of pos_bar_height to expand the top of the figure so
#' that it can capture the text of the positions.
#' @param expand_width the fraction of the width of plot to expand by so that it can
#' capture the text of the positions.
#' @param het_outline_color if this is not NA, it is the color which will
#' be used on the outside of each cell to denote whether it is heterozygous
#' or not.
#' @param plot_read_depths If this is true, this just plots a heat map
#' identical to the raster haplo that shows read depths.
#' @param annotation_columns this is a tibble that has the columns hnames,
#' column, value, width, and order.  These can be used to show different
#' characteristics in the negative x range.
#' @param annotation_rows this is a tibble that has the columns POS,
#' row, value, height, and order.  These can be used to show different
#' characteristics in the negative y range.
#' @param anno_row_start  which row below the main raster do you want the annotation rows to start?
#'
#' @export
#' @examples
#' data(RoSA)
#' data(rosa_9)
#' D <- RoSA %>%
#'      dplyr::mutate(hnames = haplo_name,
#'                    atypes = alle2)
#'
#' # make a tibble to annotate
#' pos_annot <- tibble::tibble(pos = seq(12.1e6, 12.35e6, by = 0.5e5)) %>%
#'      dplyr::mutate(pname = sprintf("%.2f", pos/1e06))
#'
#' # it turns out that these guys are ordered in RoSA how we want them to come out
#' h_ord <- unique(D$hnames)
#'
#' # color the "spring" allele yellow and the other blue
#' fcolors <- c(S = "yellow", F = "blue")
#'
#' # get positions of SNPs to highlight
#' h_pos = rosa_9 %>% dplyr::select(name, POS) %>%
#'      dplyr::mutate(POS = as.integer(POS))
#'
#' g <- haplo_raster_plot(D, h_ord, pos_annot, fcolors, h_pos)
#'
#' # if you want to express genotype posteriors by transparency
#' # you have to get those data from the unphased but imputed data
#' U <- system.file(package = "ecaRbioinf", "extdata", "greb1l-imputed-unphased.vcf.gz")
#' u <- vcfR::read.vcfR(U)
#' ut <- vcfR::vcfR2tidy(u)$gt %>%
#'      tidyr::separate(gt_GP, into = c("gp0", "gp1", "gp2"), sep = ",", convert = TRUE)
#' utj <- ut %>%
#'      dplyr::mutate(maxGP = pmax(gp0, gp1, gp2)) %>%
#'      dplyr::select(POS, Indiv, maxGP)
#'
#' D_GP <- D %>%
#'      dplyr::left_join(utj, by = c("POS", "Indiv"))
haplo_raster_plot <- function(D, h_ord = NULL, pos_annot, fcolors,
                              highlight_pos = NULL,
                              yaxis_name = "Haplotype",
                              xaxis_name = "SNP",
                              pos_bar_height = 0.10,
                              pos_bar_text_size = 1.0,
                              lo_kink_frac = 0.2,
                              hi_kink_frac = 0.2,
                              tick_frac = 0.1,
                              expand_top = 0.1,
                              expand_width = 0.0,
                              kinky_line_size = 0.1,
                              het_outline_color = NA,
                              het_outline_size = 1.0,
                              plot_read_depths = FALSE,
                              annotation_columns = NULL,
                              annotation_rows = NULL,
                              anno_row_start = 0,
                              snp_quant_tibble = NULL,
                              sq_props = list(lo = -9, hi = -1)) {

  if(!is.null(h_ord)) {
    if(any(!(h_ord %in% unique(D$hnames)))) stop("h_ord includes names not in hnames")
    D2 <- D %>%
      dplyr::filter(hnames %in% h_ord) %>%
      dplyr::mutate(Hfact = factor(hnames, levels = rev(h_ord)))
  } else {
    D2 <- D %>%
      dplyr::mutate(Hfact = factor(hnames, levels = rev(unique(D$hnames))))
  }

  D2 <- D2 %>%
    mutate(POSfact = factor(POS))

  if(!is.na(het_outline_color)) {
    D2 <- D2 %>%
      group_by(Indiv, POS) %>%
      mutate(isHet = atypes[1] != atypes[2]) %>%
      ungroup() %>%
      mutate(isHet = ifelse(isHet, "Het", "Homoz"))
  }

  # first make the base raster plot.  If we want read depths it
  # will be different
  if(plot_read_depths == FALSE) {
  base <- ggplot2::ggplot() +
    geom_raster(data = D2, aes(x = POSfact, y = Hfact, fill = atypes)) +
    scale_fill_manual(values = fcolors) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank())
  } else {
    base <- ggplot2::ggplot() +
      geom_raster(data = D2, aes(x = POSfact, y = Hfact, fill = DP)) +
      scale_fill_viridis_c() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_blank())
  }

  if(!is.na(het_outline_color)) {
    base <- base +
      geom_tile(data = D2, aes(x = POSfact, y = Hfact, colour = isHet), fill = NA, size = het_outline_size) +
      scale_colour_manual(values = c(Het = het_outline_color, Homoz = NA))
  }

  # now calculate where the position bar should go.
  topHY <- length(unique(D2$hnames)) + 0.5
  pbD <-  topHY * pos_bar_height
  pbY <-  topHY + pbD
  kinkYlo <- topHY + lo_kink_frac * pbD
  kinkYhi <- pbY - hi_kink_frac * pbD
  tickHi <- pbY + tick_frac * pbD
  topSX <- length(unique(D2$POS)) + 0.5
  hiPOS <- max(D2$POS)
  loPOS <- min(D2$POS)

  # and now make a data frame with all these parts in it
  pbtib <- pos_annot %>%
    mutate(x = topSX * (pos - loPOS) / (hiPOS - loPOS ))

  # now make a tibble with all the position information in it
  postib <- tibble::tibble(pos = sort(unique(D2$POS))) %>%
    mutate(x = topSX * (pos - loPOS) / (hiPOS - loPOS ),
           plotx = as.integer(factor(pos)))

  # now add all that stuff to the plot
  g <- base +
    geom_segment(data = pbtib, mapping = aes(x = x, xend = x, y = pbY, yend = tickHi)) +
    geom_text(data = pbtib, mapping = aes(x = x, y = tickHi, label = pname), vjust = -0.2, size = pos_bar_text_size) +
    expand_limits(y = c(0, tickHi + expand_top * pbD)) +
    expand_limits(x = c(-(expand_width * topSX), (topSX * 1 + expand_width))) +
    annotate("segment", x = min(postib$x), xend = max(postib$x), y = pbY, yend = pbY) +
    geom_segment(data = postib, mapping = aes(x = plotx, xend = plotx, y = topHY, yend = kinkYlo), size = kinky_line_size) +
    geom_segment(data = postib, mapping = aes(x = plotx, xend = x, y = kinkYlo, yend = kinkYhi), size = kinky_line_size) +
    geom_segment(data = postib, mapping = aes(x = x, xend = x, y = kinkYhi, yend = pbY), size = kinky_line_size)

  g2 <- g

  # now put some highlights around certain SNPs
  if(!is.null(highlight_pos)) {
    # first find which positions these are
    allposes <- unique(D2$POS)
    hp2 <- highlight_pos %>%
      mutate(xpos_exact_fact = factor(POS, levels = levels(D2$POSfact)),
             POSclose = unlist(lapply(POS, function(x) max(allposes[allposes <= x]))),
             xpos_close_fact = factor(POSclose, levels = levels(D2$POSfact))) %>%
      mutate(xpos_close_fact = ifelse(!is.na(xpos_exact_fact), NA, xpos_close_fact))

    # positions that are an exact match get white-and-black dashed lines around them
    # while positions that are not in the raster go between the two that they are between
    g2 <- g +
      geom_rect(data = hp2, mapping = aes(xmin = as.integer(xpos_exact_fact) - 0.5,
                                          xmax = as.integer(xpos_exact_fact) + 0.5,
                                          ymin = 0,
                                          ymax = topHY),
                fill = NA,
                size = 0.3,
                colour = "black") +
      geom_rect(data = hp2, mapping = aes(xmin = as.integer(xpos_exact_fact) - 0.5,
                                          xmax = as.integer(xpos_exact_fact) + 0.5,
                                          ymin = 0,
                                          ymax = topHY),
                fill = NA,
                size = 0.3,
                colour = "white",
                linetype = "dashed") +
      geom_segment(data = hp2, mapping = aes(x = as.integer(xpos_close_fact) + 0.5,
                                             xend = as.integer(xpos_close_fact) + 0.5,
                                             y = 0,
                                             yend = topHY),
                                             size = 0.5,
                                             colour = "darkorchid1")

  }

  # now we add the annotation columns on if we want...
  if(!is.null(annotation_columns)) {
    if(plot_read_depths == FALSE) {  # if not plotting read_depths then its all one simple layer
      Ae <- expand_anno_cols(annotation_columns, length(unique(D2$POS))) %>%
        mutate(Hfact = factor(hnames, levels = levels(D2$Hfact)))

      g3 <- g2 +
        geom_raster(data = Ae, mapping = aes(x = x, y = Hfact, fill = value))
    } else { # if plotting read depths then you have to do this hack to get discrete fills amongst your continuous fills
      Ae <- expand_anno_cols(annotation_columns, length(unique(D2$POS)))

      tmp <- tile_based_column_annotations(Ae, levels(D2$Hfact), fcolors)
      glines <- tmp$gglines
      xxx_ace_split <- tmp$xxx_ace_split

      g3 <- eval(parse(text = paste(c("g2", glines), collapse = " + ")))
    }

  } else {
    g3 <- g2
  }

  # now we add the annotation rows if we want
  if(!is.null(annotation_rows)) {
    if(plot_read_depths == FALSE) {  # if not plotting read_depths then its all one simple layer
      ARe <- expand_anno_rows(annotation_rows, length(unique(D2$hnames)), start_shift = anno_row_start) %>%
        mutate(POSfact = factor(POS, levels = levels(D2$POSfact)))

      g4 <- g3 +
        geom_raster(data = ARe, mapping = aes(x = POSfact, y = y, fill = value))
    } else {

      ARe <- expand_anno_rows(annotation_rows, length(unique(D2$hnames)), start_shift = anno_row_start)
      tmp <- tile_based_row_annotations(ARe, levels(D2$POSfact), fcolors)
      glines <- tmp$gglines
      xxx_are_split <- tmp$xxx_are_split
      g4 <- eval(parse(text = paste(c("g3", glines), collapse = " + ")))

    }

  } else {
    g4 <- g3
  }


  if(!is.null(snp_quant_tibble)) {
    sqtmp <- snp_quant_ggplot_calls(snp_quant_tibble,
                                    POS_levels = levels(D2$POSfact),
                                    top = sq_props$hi,
                                    bottom = sq_props$lo)

    xxx_SQuant <- sqtmp$xxx_SQuant

    g5 <- eval(parse(text = paste(c("g4", sqtmp$gline), collapse = " + ")))
  } else {
    g5 <- g4
  }

  g5
}



#' helper function specific to haplo_raster_plot.
#'
#' expands the annotation columns into a raster that
#' can be plotted as another layer
#' @param A a tibble with hnames, column, value, and width (which is
#' a proportion of the main raster area that you want the width of the
#' column to be).
#' @param n the width of the main raster area.  Basically the number
#' of SNPs shown.
#'
expand_anno_cols <- function(A, n) {
  tmp <- A %>%
    dplyr::mutate(cells = 1 + ceiling(n * width))

  stend <- tmp %>%
    dplyr::group_by(order) %>%
    dplyr::summarise(cells = cells[1]) %>%
    dplyr::mutate(start = dplyr::lag(cumsum(cells)),
           end = cumsum(cells)) %>%
    dplyr::mutate(start = ifelse(is.na(start), 1, start + 1)) %>%
    dplyr::select(-cells)



  dplyr::left_join(tmp, stend, by = "order") %>%
    dplyr::group_by(column) %>%
    dplyr::do(tidyr::expand(., tidyr::nesting(hnames, column, value), x = .$start[1]:.$end[1])) %>%
    dplyr::mutate(x = -x)
}







#' helper function specific to haplo_raster_plot.
#'
#' expands the annotation rows into a raster that
#' can be plotted as another layer
#' @param A a tibble with POS, row, value, and height (which is
#' a proportion of the main raster-area height that you want the height of the
#' row to be).
#' @param n the height of the main raster area.  Basically the number
#' of hnames shown.
#' @param start_shift How many rows down do we want to start these guys?
#'
expand_anno_rows <- function(A, n, start_shift = 0) {
  tmp <- A %>%
    dplyr::mutate(cells = 1 + ceiling(n * height))

  stend <- tmp %>%
    dplyr::group_by(order) %>%
    dplyr::summarise(cells = cells[1]) %>%
    dplyr::mutate(start = dplyr::lag(cumsum(cells)),
                  end = cumsum(cells)) %>%
    dplyr::mutate(start = ifelse(is.na(start), 1, start + 1)) %>%
    dplyr::select(-cells) %>%
    dplyr::mutate(start = start + start_shift,
                  end = end + start_shift )



  dplyr::left_join(tmp, stend, by = "order") %>%
    dplyr::group_by(row) %>%
    dplyr::do(tidyr::expand(., tidyr::nesting(POS, row, value), y = .$start[1]:.$end[1])) %>%
    dplyr::mutate(y = -y)
}




#' This function creates layers using geom_tile from *expanded* column_annotations
#'
#' If fill is mapped to a continuous variable, but you still want your annotations
#' to have discrete values, this is how you have to do it.
#'
#' This thing returns a vector of ggplot commands to add the layers
#' @param ace the expanded column annotations data frame.  It must have
#' the columns: hnames, column, value, x.
#' @param h_ord how the hnames in D are supposed to be ordered
#' @param fcolors a named vector of colors to use for the discrete values
tile_based_column_annotations <- function(ace, h_ord, fcolors) {

  # get y-values according to how hnames is ordered, and add a color column on there according to column and value
  ace2 <- ace %>%
    dplyr::mutate(y = as.integer(factor(hnames, levels = h_ord))) %>%
    dplyr::mutate(color = fcolors[value])


  # then split it on column and value
  ace_split <- split(ace2, list(ace2$column, ace2$value))

  # and then remove all those combinations of column and value that have no rows...
  ace_split <- ace_split[unlist(lapply(ace_split, nrow)) > 0]


  # then make a bunch of geom_tile lines that we will want to bung together with pluses and eval:
  gglines <- unlist(
    lapply(seq_along(ace_split), function(n) {
      paste0("ggplot2::geom_tile(data = xxx_ace_split[[", n, "]], mapping = ggplot2::aes(x = x, y = y, width = 1, height = 1), fill = \"", ace_split[[n]]$color[1] , "\")")
    })
  )

  list(gglines = gglines, xxx_ace_split = ace_split)

}




#' This function creates layers using geom_tile from *expanded* row_annotations
#'
#' If fill is mapped to a continuous variable, but you still want your annotations
#' to have discrete values, this is how you have to do it.
#'
#' This thing returns a vector of ggplot commands to add the layers
#'
#' This assumes that the POS values in are are just the ones that should be there.
#'
#' @param are the expanded row annotations data frame.  It must have
#' the columns: POS, row, value, y.
#' @param pos_values the levels of the POS values to make them a factor from
#' @param fcolors a named vector of colors to use for the discrete values
tile_based_row_annotations <- function(are, pos_values, fcolors) {

  # get x-values according to the POS, and add a color column on there according to column and value
  are2 <- are %>%
    dplyr::mutate(x = as.integer(factor(POS, levels = pos_values))) %>%
    dplyr::mutate(color = fcolors[value])


  # then split it on column and value
  are_split <- split(are2, list(are2$row, are2$value))

  # and then remove all those combinations of column and value that have no rows...
  are_split <- are_split[unlist(lapply(are_split, nrow)) > 0]


  # then make a bunch of geom_tile lines that we will want to bung together with pluses and eval:
  gglines <- unlist(
    lapply(seq_along(are_split), function(n) {
      paste0("ggplot2::geom_tile(data = xxx_are_split[[", n, "]], mapping = ggplot2::aes(x = x, y = y, width = 1, height = 1), fill = \"", are_split[[n]]$color[1] , "\")")
    })
  )

  list(gglines = gglines, xxx_are_split = are_split)

}



#' this is a function to plot quantitative information about each SNP
#'
#' This bad boy will plot a background tile that covers anything.  It is rimmed with
#' a line of "rim_color" and it is filled with "background_color" and then it throws down
#' the SNP-specific information as a height bar of height "value" (from the tibble)
#'
#' This thing returns a vector of ggplot calls that can be run through eval, and also
#' the tibble to stuff into it.  (Doing it
#' this way lets us mix continuous and discrete color scales....)
#'
#' @param Q the tibble that has POS and value.
#' @param POS_levels the levels of the POSes in the current raster plot.  This will
#' determine the x position when we coerce POS to a factor.
#' @param top the row number of the top of the background (corresponds to value = 1)
#' @param bottom the row number of the bottom of the background (corresponds to value = 0)
#' @param background_color string which is the color desired for the background
#' @param color string which is the color desired for the value bar
#' @param rim_color color for the line around the box of information
snp_quant_ggplot_calls <-function(Q,
                                  POS_levels,
                                  top, bottom, background_color = "white", color = "black", rim_color = "black") {

  Q2 <- Q %>%
    mutate(POSfact = as.integer(factor(POS, levels = POS_levels))) %>%
    mutate(ylo = bottom,
           yhi = bottom + (top - bottom) * value)


  bg <- paste0("annotate(\"rect\", xmin = 0, xmax = ",length(POS_levels), ", ymin = ", bottom, ", ymax = ", top, ", colour = \"", rim_color, "\", fill = \"", background_color, "\", size = 0.1)")

  vals <- paste0("geom_rect(data = xxx_SQuant, mapping = aes(xmin = POSfact - 0.5, xmax = POSfact + 0.5, ymin = ylo, ymax = yhi), fill = \"", color, "\", colour = NA)")

  list(gline = paste(bg, vals, sep = " + "),
       xxx_SQuant = Q2)
  }
