#' match locations on a short query sequence to the genomic coordinates of a reference
#'
#' This is particularly useful when you mapped to an early version of a genome, and
#' then designed assays for SNPs at positions X, Y, and Z in a, say, 500 bp
#' fragment.  But now you want to know the coordinates of X, Y, and Z relative
#' to an updated version of the genome.  So, you map the 500 bp to the updated
#' version of the genome, only to find that there is a long and complicated
#' CIGAR string that links the two.  This parses that CIGAR string to give you
#' the genomic coordinates of X, Y, and Z relative to the new, updated reference
#' genome.
#'
#' This is definitely not written for speed.  It would not be hard to write something
#' super fast in a few lines of C++, but I wanted to do something that is easy
#' to look at.
#'
#' For the most part, this can be done using cumulative sums on whether or not
#' each different CIGAR operator consumes the query or the reference.
#'
#' It looks like the same or a similar thing can be done in the BioConductor
#' package function,
#' GenomicAlignments::mapToAlignments, but I find that stuff pretty opaque
#' (so many different S4 classes....I suppose if one was going to take a quarter
#' to learn all about it, it would be good for doing bioinformatics.  But, for
#' an ecologist, it is, I think, way better to spend that quarter using how to use
#' the tidyverse well.)
#' @param cigar the cigar string
#' @param start the base-1 left end mapping position of the query on the reference
#' @param rel_pos the base-1 positions (1 is the start of the query sequence) whose
#' corresponding positions in the reference are desired.  If NULL, all positions in the
#' query are returned.  It is very important to note that the position here is
#' relative to the original query sequence, upon which hard-clipping has not been
#' done, for that reason, hard-clipping consumes this original query sequence.
#' @param strand  This should be "+" or "-".  "+" indicates an alignment to the
#' forward strand and it means that the rel_pos on the query is correct as it
#' stands.  If strand is "-" then the function takes care of that by computing
#' a new relative position that starts counting from the other end of the query
#' strand.  That is why query length is required.
#' @param query_length The length of the original query sequence.
#' @export
genome_coords_from_query_and_cigar <- function(cigar, start, rel_pos = NULL, strand = NULL, query_length = NULL) {
  # first, parse the cigar string into a tibble and add the incrementers on it
  cig <- tibble::tibble(
    len = as.integer(head(str_split(cigar, "[MIDNSHP=X]")[[1]], -1)),
    op = str_split(cigar, "[0-9]+")[[1]][-1]
  ) %>%
    dplyr::left_join(cigar_increments_orig_query, by = "op")

  # then replicate all the lines as indicated and compute the query position and the
  # reference position for each row
  n <- nrow(cig)
  cigr <- cig[rep(1:n, cig$len), ] %>%
    mutate(query_pos = cumsum(query_increment),
           reference_pos = start - 1 + cumsum(reference_increment))

  ret <- cigr

  # so, now that big table is made, we just need to pick out the correct relative positions from it
  if(!is.null(rel_pos)) {
    if(!is.null(strand)) {
      if(strand == "-") {
        if(is.null(query_length)) stop("You must supply query_length of strand is \"-\"")
        rel_pos <- quuery_length -  rel_pos + 1
      }
    }
    ret <- cigr %>%
      dplyr::filter(query_pos %in% rel_pos)
  }

  ret
}
