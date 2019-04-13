
#' convert phased variants in VCF to tidy data frame of haplotypes
#'
#' If the ancestral state is listed this gets put out there too
#' @param V the path to the VCF.  Note that the allele separator has to be "|" for
#' everything, as it has to be phased.
#' @param Anc if NULL (the default) this is ignored.  If non-null, this is the name
#' to be given to the ancestral sequence.  If used, this requires that the VCF has
#' an AA field within the INFO field.
#' @return A list with at least two components: \code{fix}: the fixed infor about the variants
#' from the VCF, and \code{tidy} a tidy data frame with columns: ChromKey, POS, Indiv,
#' haplo (a or b), allele.
#' If Anc is given then that will be a sample with an a-haplo in \code{tidy}.  However
#' If Anc is given, then the return list will have a third component named
#' avd, which is just
#' like tidy, except that it includes an extra column \code{anc_vs_derived} that
#' has an A if the allele is ancestral (i.e. same as the ancestral sequence), a D
#' if it is derived (different from the ancestral sequence).
#' (And also, this data frame does not include the ancestral sequence as an Indiv.)
#' @export
#' @examples
#' V <- system.file(package = "ecaRbioinf", "extdata", "greb1l-imp-phased-with-anc.vcf.gz")
#' haps <- vcf_haplos2tidy(V, Anc = "Coho")
vcf_haplos2tidy <- function(V, Anc = NULL) {
  # read and tidy VCF
  v <- vcfR::read.vcfR(V)
  vt <- vcfR::vcfR2tidy(v)

  # get phased haplotypes a and b
  vtp <- vt$gt %>%
    tidyr::separate(gt_GT_alleles, into = c("a", "b"), sep = "\\|") %>%
    dplyr::select(ChromKey, POS, Indiv, a, b) %>%
    tidyr::gather(key = "haplo", value = "allele", a, b)

  # add ancestral sequence if requested
  if(!is.null(Anc)) {
    Anc <- as.character(Anc)
    stopifnot(length(Anc) == 1)

    if(!("AA" %in% names(vt$fix))) {
      stop("It appears your VCF does not have an AA field in INFO while you are reqeusting an ancestral sequence.")
    }

    vtpAA <- vt$fix %>%
      dplyr::mutate(Indiv = Anc, haplo = "a", allele = AA) %>%
      dplyr::select(ChromKey, POS, Indiv, haplo, allele) %>%
      dplyr::mutate(allele = ifelse(allele == "N", NA, allele))  # if no allele (i.e. "N") is available from the ancestral sequence, then call it NA

    vtp <- dplyr::bind_rows(vtp, vtpAA)
  }

  # arrange appropriately
  vtp <- vtp  %>%
    dplyr::arrange(Indiv, haplo, ChromKey, POS)

  ret <- list(fix = vt$fix, tidy = vtp)

  # now, deal with classifying alleles as A, D, or M
  if(!is.null(Anc)) {
    # get the ancestral alleles and then left_join them onto everyone
    ancy <- vtp %>%
      dplyr::filter(Indiv == Anc) %>%
      dplyr::select(ChromKey, POS, allele) %>%
      dplyr::rename(anc_allele = allele)

    vta <- dplyr::left_join(vtp, ancy, by = c("ChromKey", "POS")) %>%
      dplyr::mutate(anc_vs_derived = ifelse(allele == anc_allele, "A", "D")) %>%
      dplyr::select(-allele, -anc_allele) %>%
      dplyr::filter(Indiv != Anc)


    ret <- list(fix = vt$fix, tidy = vtp, avd = vta)

  }

  ret
}
