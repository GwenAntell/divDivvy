#' Find unique (taxon) occurrence records
#'
#' Subset a dataset to unique spatial localities or locality-taxon combinations.
#'
#' @inheritParams sdsumry
#' @param siteId The name or numeric position of the column containing
#' identifiers for unique spatial sites, e.g. raster cell names.
#' If `siteId` supplied, `xy` is ignored. One of `siteId` or `xy` is required.
#' @param na.rm Should records missing spatial information be removed?
#' Default is yes.
#'
#' @return An object with the same class and columns as `dat`, containing the
#' subset of rows representing unique localities (if `siteId` supplied),
#' coordinates (if `xy` supplied), or taxon-site combinations
#' (if `taxVar` supplied). The first record at each spatial locality is retained,
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
#' # rows 4 and 5 removed as duplicates of 1 and 2, respectively
#' obs
#' uniqify(obs, taxVar = 3, xy = 1:2)
#'
#' # using taxon identifications or other third variable is optional
#' uniqify(obs, xy = c('x', 'y'))
#'
#' # caution - data outside the taxon and occurrence variables
#' # will be lost where associated with duplicate occurrences
#' obs$notes <- letters[11:20]
#' uniqify(obs, taxVar = 3, xy = 1:2)
#' # the notes 'n' and 'o' are absent in the output data
#'
#' @export

uniqify <- function(dat, taxVar = NULL, xy = NULL, siteId = NULL,
                    na.rm = TRUE){
  if (is.null(siteId) & is.null(xy)) stop('supply column name(s) for spatial location')

  # if siteId is supplied, always use that over coordinates
  # to minimise problems finding identical floating-point numbers
  if (! is.null(siteId)){
    args <- siteId
  } else {
    if (! is.null(xy)){
      args <- xy
    }
  }
  if (! is.null(taxVar)){
    args <- c(taxVar, args)
  }

  if (na.rm){
    compl <- stats::complete.cases(dat[,args])
    dat <- dat[compl, ]
  }
  dupes <- duplicated( dat[,args] )
  dat[ !dupes, ]
}
