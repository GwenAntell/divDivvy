#' Find unique (taxon) occurrence records
#'
#' Subset a dataset to unique spatial localities or locality-taxon combinations.
#'
#' @inheritParams sdsumry
#' @param na.rm Should records missing spatial information be removed?
#' Default is yes.
#'
#' @return An object with the same class and columns as `dat`, containing the
#' subset of rows representing unique coordinates (if only `xy` supplied)
#' or unique taxon-site combinations (if `taxVar` is also supplied).
#' The first record at each spatial locality is retained,
#' or if `taxVar` is specified, the first record of each taxon at a locality.
#'
#' @examples
#' # generate occurrence data
#' x  <- rep(1, 10)
#' y  <- c(rep(1, 5), 2:6)
#' sp <- c(rep(letters[1:3], 2),
#'         rep(letters[4:5], 2))
#' obs <- data.frame(x, y, sp)
#'
#' # compare original and unique datasets:
#' # rows 4 and 5 removed as duplicates of rows 1 and 2, respectively
#' obs
#' uniqify(obs, taxVar = 3, xy = 1:2)
#'
#' # using taxon identifications or other third variable is optional
#' uniqify(obs, xy = c('x', 'y'))
#'
#' # caution - data outside the taxon and occurrence variables
#' # will be lost where associated with duplicate occurrences
#' obs$notes <- letters[11:20]
#' uniqify(obs, 1:2, 3)
#' # the notes 'n' and 'o' are absent in the output data
#'
#' @export

uniqify <- function(dat, xy, taxVar = NULL, na.rm = TRUE){
  args <- xy
  if (! is.null(taxVar)){
    args <- c(taxVar, args)
  }
  if (na.rm){
    compl <- stats::complete.cases( dat[,args] )
    dat <- dat[compl, ]
  }
  dupes <- duplicated( dat[,args] )
  dat[ !dupes, ]
}
