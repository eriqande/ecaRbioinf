#' Find and return the substitutions between two parallel fasta files (each with only a single sequence)
#'
#' This is for a case like this: you have a fasta file for one chromosome in a species,
#' for example Chinook salmon.  You also have another fasta file which holds the bases
#' in coho salmon at the homologous positions in Chinook salmon.  (In other words, the
#' coho fasta is exactly the same length as the Chinook fasta.)  You want to find all the
#' positions where you have a different base in Chinook than you do in coho, and you want
#' to return those all in a nice data frame.
#'
#' Note that upper and lower case letters are treated as the same.  So, for example
#' ref = A and anc = a will not be counted as a substitution.  Output will be always
#' uppercase.
#'
#' Additionally, if either sequence has an N, that will not be counted as a difference
#' or a substitution.  In fact, all bases that are not A,C,G, or T will be discarded
#' and not recorded as a difference.  (In other words, we are just going for substitutions
#' here).
#'
#'
#' Note, the chromosome name (CHROM in the output) will be taken from the reference
#' fasta (ref).
#'
#' Note, if the sequence name in ref is of the form >chrom:p1-p2, where p1 and p2
#' are integers, then the positions
#' of all the reported differences will be shifted to start at p1.
#'
#' It was surprising to me that these operations could be done so quickly and easily
#' in R.
#'
#' @param ref path to a fasta (can be gzipped) that contains the "reference" sequence.
#' @param anc path to a fasta (can be gzipped) which represents that homologous sequence,
#' which will often be the "ancstral" sequence which has been aligned to the exact positions
#' in the reference sequence.
#' @param OutputBigTable if TRUE, then this will return a list, element "ret" of which is
#' the standard return value, and element "raw_tib" of which is the full tibble with the
#' full sequences.
#' @export
#' @examples
#' ref <- system.file(package = "ecaRbioinf", "extdata", "NC_037124.1_fragment_chinook.fna.gz")
#' anc <- system.file(package = "ecaRbioinf", "extdata", "NC_037124.1_fragment_coho.fna.gz")
#' fdd <- fasta_dna_diffs(ref, anc)
fasta_dna_diffs <- function(ref, anc, OutputBigTable = FALSE) {

    rs <- scan(ref, what = "character")
    as <- scan(anc, what = "character")

    # check to make sure there is only a single sequence in each:
    r_ncounts <- sum(stringr::str_detect(rs, "^>"))
    a_ncounts <- sum(stringr::str_detect(as, "^>"))

    if(r_ncounts > 1) {
      stop("Apparently multiple sequences in ref")
    }
    if(a_ncounts > 1) {
      stop("Apparently multiple sequences in anc")
    }

    # check to make sure that the first element is the chromsome name
    if(!(stringr::str_detect(rs[1], "^>")[1])) {
      stop("Sequence name in ref doesn't start with \">\" like it should in a fasta file.")
    }
    if(!(stringr::str_detect(as[1], "^>")[1])) {
      stop("Sequence name in ref doesn't start with \">\" like it should in a fasta file.")
    }

    # parse the chromosome name in ref in order to infer the starting position
    start_pos <- 1L  # this is the default
    tmp <- stringr::str_split(rs[1], ":")[[1]]
    chrom_name <- stringr::str_replace(tmp[1], ">", "")
    if(!is.na(tmp[2])) {
      start_pos <- as.integer(stringr::str_split(tmp[2], "-")[[1]][1])
      if(is.na(start_pos)) {
        stop("Sequence name had a colon and I expected a parseable integer start location.  Something is wrong. Offending sequence name is ", rs[1])
      }
    }
    message("Sequence starting position: ", start_pos)


    # now, we turn the sequences into character vectors.  One base for each element:
    rseq <- toupper(strsplit(paste(rs[-1], collapse = ""), "")[[1]])
    aseq <- toupper(strsplit(paste(as[-1], collapse = ""), "")[[1]])

    if(length(rseq) != length(aseq)) {
      stop("Sequence length of ref and anc are not the same!")
    }

    # now we make a tibble of these.
    st <- tibble::tibble(CHROM = chrom_name,
                         POS = seq(start_pos, length.out = length(rseq)),
                         f_ref = rseq,
                         f_anc = aseq)


    # now we filter things down to only the substitutions
    bases <- c("A", "C", "G", "T")
    tmp <- ret <- st %>%
      dplyr::filter( (f_ref %in% bases) & (f_anc %in% bases))

    subs <- tmp %>%
      dplyr::filter(f_ref != f_anc)

    message("Percentage substitutions: ", nrow(subs)/nrow(tmp))

    if(OutputBigTable == TRUE) {
      ret <- list(ret = subs, raw_tib = st)
    } else {
      ret <- subs
    }

    ret
}
