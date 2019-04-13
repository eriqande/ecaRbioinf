#' insert variants into sequences
#'
#' This assumes a fasta file that is available that includes
#' all the sequence that you want to propagate the variants
#' into.
#' @param fasta path to the fasta file.  Currently it has to have only a
#' single segment in it, and you should already have extracted the
#' part of the sequence that you want to propagate the variants into.
#' This file can be gzipped. The name of the segment has to
#' have the start position in it (relative to the start of the chromosome
#' as presumably used in the POS in V)
#' @param variant_extent_overhang If this is NULL then the full sequence from
#' fasta is used.  Otherwise, make it an integer, in which case that is the
#' number of base pairs from min and max of V$POS that are included in the sequence.
#' @param V a tidy data frame like that produced by vcf_haplos2tidy()
#' @param lo 1-based genomic coordinate to start from, if desired.
#' @param hi 1-based genomic coordinate to end at. Note that lo and hi must be contained
#' within the fasta.
#' @export
#' @examples
#' data(haps)
#' V <- haps$tidy
#' fasta <- system.file("extdata", "NC_037124.1_fragment_chinook.fna.gz", package = "ecaRbioinf")
#' lo <- 12260034
#' hi <- 12289484
insert_variants_into_sequences <- function(V, fasta, file = "tmp.txt",
                                           lo = NA, hi = NA) {

  rs <- scan(fasta, what = "character")

  # check to make sure there is only a single sequence in each:
  r_ncounts <- sum(stringr::str_detect(rs, "^>"))
  if(r_ncounts > 1) {
    stop("Apparently multiple sequences in ref")
  }

  # check to make sure that the first element is the chromsome name
  if(!(stringr::str_detect(rs[1], "^>")[1])) {
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

  # break it into a vector
  rseq <- toupper(strsplit(paste(rs[-1], collapse = ""), "")[[1]])


  # now some checking
  fhi_pos <- start_pos + length(rseq) - 1

  # now trim the sequence down to the desired lo and hi with some error checking
  if(is.na(lo)) lo <- start_pos
  if(is.na(hi)) hi <- fhi_pos

  if(lo < start_pos || hi > fhi_pos) stop("lo or hi out of range of sequence")

  rseq <- rseq[-c(
    0:(lo - start_pos),
    (hi - start_pos + 2):(fhi_pos - start_pos + 2)
  )]
  start_pos <- lo
  fhi_pos <- hi


  if(min(V$POS) < start_pos) {
    message("fasta file start point ",
            start_pos,
            " is larger than the lowest POS in V: ",
            min(V$POS),
            ".  Discarding POSes below ",
            start_pos)
  }
  if(max(V$POS) > fhi_pos) {
    message("fasta file end point ",
            fhi_pos,
            " is smaller than the largest POS in V: ",
            max(V$POS),
            ".  Discarding POSes above ",
            fhi_pos)
  }

  # now, filter the POSes and mutate in a new position for each marker and arrange appropriately
  Vsplit <- V %>%
    filter(POS >= start_pos, POS <= fhi_pos) %>%
    mutate(newpos = POS - start_pos + 1) %>%
    arrange(Indiv, haplo, newpos) %>%
    split(., f = list(.$Indiv, .$haplo))



  # note that this is trivial to parallelize with mclapply if desired
  inserted_seqs <- lapply(Vsplit, function(x) {
    ret <- rseq
    ret[x$newpos] <- x$allele
    ret
  })


  # then write that to a file
  if(file.exists(file)) file.remove(file)

  dump <- lapply(names(inserted_seqs), function(n) {
    cat(">", n, "\n", file = file, append = TRUE, sep = "")
    cat(inserted_seqs[[n]], "\n", file = file, append = TRUE, sep = "")
  })

}
