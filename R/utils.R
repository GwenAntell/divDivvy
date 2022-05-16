#' Find unique (taxon) occurrence records
#'
#' Subset a dataset to unique spatial localities or locality-taxon combinations.
#'
#' @inheritParams sdsumry
#' @param siteId The name or numeric position of the column containing
#' identifiers for unique spatial sites, e.g. raster cell names.
#' If `siteId` supplied, `xy` is ignored.
#' @param na.rm Should records missing spatial information be removed?
#' Default is yes.
#'
#' @return An object with the same class and columns as `dat`, containing the
#' subset of rows representing unique localities (if `siteId` supplied),
#' coordinates (if `xy` supplied), or taxon-site combinations
#' (if `taxVar supplied`). The first record at each spatial locality is retained,
#' or if `taxVar` is specified, the first record of each taxon at a locality.
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
