

#' Add dna differences from fasta files to phased tidy vcfR data frame
#'
#' Here is the situation where you might use this: you have phased sequence
#' data from chinook salmon, and you also have the ancestral state from coho
#' for each of those sites.  But now you want to add into that all the sites
#' where there are differences between coho and Chinook, although Chinook has
#' no variants at these sites.
#'
#' Note that a single CHROM is allowed and it is assumed to be the same for HT and ref
#' and anc.  No checking of this is done at the moment.
#' @param HT tidy data of haplotypes, like that coming from vcf_haplos2tidy().
#' @param Anc the name of the ancestral sequence in HT
#' @inheritParams fasta_dna_diffs
#' @export
#' @examples
#' # get the VCF first
#' V <- system.file(package = "ecaRbioinf", "extdata", "greb1l-imp-phased-with-anc.vcf.gz")
#' haps <- vcf_haplos2tidy(V, Anc = "Coho")
#'
#' # whittle it down to the RoSA region
#' haps$fix <- haps$fix %>% dplyr::filter(POS > 12.05e6, POS < 12.4e6)
#' haps$tidy <- haps$tidy %>% dplyr::filter(POS > 12.05e6, POS < 12.4e6)
#'
#' # get the fasta
#' ref <- system.file(package = "ecaRbioinf", "extdata", "NC_037124.1_fragment_chinook.fna.gz")
#' anc <- system.file(package = "ecaRbioinf", "extdata", "NC_037124.1_fragment_coho.fna.gz")
#'
#' # add the substitutions
#'
#' combo <- add_substitutions(haps, "Coho", ref, anc)
add_substitutions <- function(HT, Anc, ref, anc) {

  # get the DNA diffs and filter out any that are outside of the POS range of HT
  fdd <- fasta_dna_diffs(ref, anc)

  fd2 <- fdd %>%
    dplyr::filter(POS > min(HT$tidy$POS), POS < max(HT$tidy$POS))

  # remove any substitutions that coincide with variation
  fd3 <- fd2 %>%
    dplyr::anti_join(HT$tidy, by = "POS")

  expos <- expand.grid(POS = fd3$POS, Indiv = unique(HT$tidy$Indiv), haplo = unique(HT$tidy$haplo)) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(POS = as.integer(as.character(POS)),
                  Indiv = as.character(Indiv),
                  haplo = as.character(haplo))

  combined <- expos %>%
    dplyr::left_join(fd3, by = "POS") %>%
    dplyr::mutate(allele = ifelse(Indiv == Anc, f_anc, f_ref)) %>%
    dplyr::mutate(ChromKey = 1) %>%
    dplyr::select(-f_ref, -f_anc, -CHROM) %>%
    dplyr::select(ChromKey, dplyr::everything())

  # finally bind those together and re-order
  ret <- dplyr::bind_rows(HT$tidy, combined) %>%
    dplyr::arrange(Indiv, haplo, POS)

  ret
}
