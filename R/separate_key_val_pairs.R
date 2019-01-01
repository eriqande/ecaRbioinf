

#' separate a column of multiple key=value pairs into separate columns headed by the keys
#'
#' This is an ugly problem when we have different keys in different rows.  Here I am
#' going to try an easy-to-code solution that not be super efficient, but interesting
#' to see, nonetheless.
#'
#' The slow part is actually just parsing things out with the regexes.  Beyond that
#' it is all pretty fast.
#'
#' Note that this does not currently convert columns to different variable types afterward.
#' That is up to the user.
#'
#' @param D a data frame
#' @param column a string-name of the column to separate
#' @param key_regex reg-ex to pick out they keys.  Default is for =/;
#' @param val_regex reg-ex to pick out the vals. Default is for =/;
#' @param remove if TRUE the column is removed
#' @export
#' @examples
#' data(snp_anno_intersect)
#' D <- snp_anno_intersect
#' D2 <- separate_key_value_pairs(D, "biggie")
separate_key_value_pairs <- function(D,
                                     column,
                                     key_regex = "([a-zA-Z_0-9]+)=",
                                     val_regex = "=([^;]*)",
                                     remove = TRUE
) {


  keys_list <- stringr::str_match_all(D[[column]], key_regex)
  vals_list <- stringr::str_match_all(D[[column]], val_regex)

  list_of_named_vals <- lapply(seq_along(keys_list), function(i) {
    x <- vals_list[[i]][, 2];
    names(x) <- keys_list[[i]][,2];
    x
    })

  all_names <- unique(names(unlist(list_of_named_vals)))

  # now lapply over those and  add values into those
  full_list <- lapply(list_of_named_vals, function(x) {
    y <- rep(NA_character_, length(all_names))
    names(y) <- all_names
    y[names(x)] <- x
    y
  })

  # now, make a matrix of it and then turn it into a tibble
  mat <- matrix(unlist(full_list), byrow = TRUE, ncol = length(all_names))

  new_columns <- tibble::as_tibble(mat) %>%
    setNames(all_names)

  # check for name overlap between new columns and old
  lappy <- names(new_columns) %in% names(D)

  while(any(lappy)) {
    warning("Keys of new columns overlap with old. Will be pre-pended with \"x.\": ",
            paste(names(new_columns)[lappy], collapse = ", "))
    names(new_columns)[lappy] <- paste0("x.", names(new_columns)[lappy])
    lappy <- names(new_columns) %in% names(D)
  }

  # then bind_cols them back to the original
  if(remove == TRUE) {
    ret <- D %>%
      dplyr::select(- !! rlang::sym(column)) %>%
      dplyr::bind_cols(new_columns)
  } else {
    ret <- dplyr::bind_cols(D, new_columns)
  }

  ret
}



# this is a SLOW version
separate_key_value_pairs_SLOW <- function(D,
                                     column,
                                     key_regex = "([a-zA-Z_0-9]+)=",
                                     val_regex = "=([^;]*)",
                                     remove = TRUE
) {
  Dlist <- split(D, f = 1:nrow(D))

  new_columns <- lapply(Dlist, function(x) {
    keys <- stringr::str_match_all(x[column], key_regex)[[1]][,2]
    vals <- stringr::str_match_all(x[column], val_regex)[[1]][,2]
    tmp <- matrix(vals, nrow = 1)
    colnames(tmp) <- keys
    tibble::as_tibble(tmp)
  }) %>%
    dplyr::bind_rows()

  # check for name overlap between new columns and old
  lappy <- names(new_columns) %in% names(D)

  while(any(lappy)) {
    warning("Keys of new columns overlap with old. Will be pre-pended with \"x.\": ",
            paste(names(new_columns)[lappy], collapse = ", "))
    names(new_columns)[lappy] <- paste0("x.", names(new_columns)[lappy])
    lappy <- names(new_columns) %in% names(D)
  }

  # then bind_cols them back to the original
  if(remove == TRUE) {
    ret <- D %>%
      dplyr::select(- !! sym(column)) %>%
      dplyr::bind_cols(new_columns)
  } else {
    ret <- dplyr::bind_cols(D, new_columns)
  }

  ret
}


