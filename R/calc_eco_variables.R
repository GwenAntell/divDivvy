
rangeSizer <- function(coords){
  latDiff <- max(coords[,2]) - min(coords[,2])
  latRng <- abs(latDiff)
  pts <- sf::st_as_sf(coords, coords = 1:2, crs = 'epsg:4326')
  ptsGrp <- sf::st_union(pts)
  cntr <- unlist( sf::st_centroid(ptsGrp) )
  gcdists <- sf::st_distance(pts) # returns units-class object (m)
  gcdists <- units::set_units(gcdists, 'km')
  gcMax <- max(gcdists)
  mst <- vegan::spantree( units::drop_units(gcdists) )
  agg <- sum(mst$dist)
  out <- cbind('centroidLng' = cntr[1],
               'centroidLat' = cntr[2],
               'latRange' = latRng,
               'greatCircDist' = gcMax,
               'minSpanTree' = agg
  )
  return(out)
}

#' Calculate basic spatial coverage and diversity metrics
#'
#' Summarise the geographic scope and position of occurrence data, and
#' optionally estimate richness and evenness
#'
#' \code{sampMeta} compiles a variety of metadata about a (sub)sample,
#' before or after spatial subsampling. The function counts the number
#' of unique spatial sites, collections (if requested), and taxa, and
#' calculates the centroid coordinates, latitudinal range (degrees),
#' great circle distance (km), and summed minimum spanning tree length (km)
#' for occurrences.
#'
#' If \code{quotaQ} is supplied, \code{sampMeta} rarefies richness at the
#' given coverage value and returns the point estimate of richness (Hill number 0)
#' and its 95% confidence interval as well as estimates of evenness (Pielou's J)
#' and sample coverage (Good's *u*). If \code{quotaN} is supplied,
#' \code{sampMeta} rarefies richness to the given number of occurrence counts
#' and returns the point estimate of richness and its 95% confidence interval.
#'
#' Coverage-based and classical rarefaction are both calculated with
#' \code{iNEXT::estimateD} internally. For details, such as how diversity
#' is extrapolated if a sample is insufficient to achieve a specified
#' rarefaction level, consult Chao and Jost (2012) and Hsieh *et al.* (2016).
#'
#' @param dat A \code{data.frame} or \code{matrix} containing taxon names,
#' coordinates, and any associated variables.
#' @param taxVar The name or numeric position of the column containing
#' taxonomic identifications.
#' @param xy A vector of two elements, specifying the name or numeric position
#' of the columns containing longitude and latitude coordinates, respectively.
#' @param collections The name or numeric position of the column containing
#' unique collection IDs, e.g. 'collection_no' in PBDB data downloads.
#' @param quotaQ A numeric value for the coverage (quorum) level at which to
#' perform coverage-based rarefaction (shareholder quorum subsampling).
#' @param quotaN A numeric value for the quota of taxon occurrences to subsample
#' in classical rarefaction.
#' @param omitDom If \code{omitDom = TRUE} and \code{quotaQ} or \code{quotaN}
#' is supplied, remove the most common taxon prior to rarefaction.
#'
#' @export
#'
#' @references
#'
#' \insertRef{Chao2012}{divvy}
#'
#' \insertRef{Hsieh2016}{divvy}

sampMeta <- function(dat, taxVar, xy,
                     collections = NULL,
                     quotaQ = NULL, quotaN = NULL,
                     omitDom = FALSE){
  out <- c()
  # n. collections (PBDB data) may be useful to know if rarefying
  if ( !is.null(collections) ){
    collSamp <- unique( dat[,collections] )
    out <- cbind(out, 'nColl' = length(collSamp))
  }
  # comb out any duplicate occurrences of a taxon w/in single site.
  # do this after counting collections in case a single taxon
  # at a single site is recorded in multiple collections
  dat <- unique( dat[,c(taxVar, xy)] )

  dat[,'site'] <- paste(dat[, xy[1] ], dat[, xy[2] ], sep = '/')
  siteIds <- unique( dat[,'site'] )
  nSite <- length(siteIds)
  out <- cbind(out, 'nSite'= nSite)

  # run range size fcn as if all occs were from a single taxon
  # duplicate localities affect only occ count, nothing else
  sites <- unique( dat[,xy] )
  spatMeta <- rangeSizer(coords = sites)
  out <- cbind(out, spatMeta)

  # diversity metrics; optional coverage and/or classical rarefaction
  taxSamp <-  unique( dat[,taxVar] )
  s <- length(taxSamp)
  out <- cbind(out, 'nTax' = s) # raw count of taxa
  freqs <- table( dat[,taxVar] )
  freqOrdr <- sort( as.numeric(freqs), decreasing = TRUE)

  # Pielou's J evenness metric (calculate before omitting dom)
  pTax <- freqOrdr / sum(freqOrdr)
  h <- -sum(pTax * log(pTax))
  # vegan::diversity('shannon', base = exp(1))
  j <- h / log(s)

  if (omitDom == TRUE){
    dom <- names( which.max(freqs) )
    domRows <- dat[,taxVar] == dom
    dat <- dat[ !domRows, ]
  }
  if ( !is.null(quotaQ) ){
    sqsFull <- iNEXT::estimateD(list( c(nSite, freqOrdr) ),
                                datatype = 'incidence_freq',
                                base = 'coverage', level = quotaQ # , conf=NULL # 95% CI or none
                                )
    # q0 = richness, q1 = exp of Shannon's entropy index,
    # q2 = inverse of Simpson's concentration index
    sqsRich <- sqsFull[sqsFull$order == 0, c('SC','qD','qD.LCL','qD.UCL')]
    names(sqsRich) <- c('u','SQSdiv','SQSlow95','SQSupr95')
    # return sample coverage, richness estimate and 95% CI
    out <- cbind(out, 'evenness' = j, sqsRich)
  }
  if ( !is.null(quotaN) ){
    crFull <- iNEXT::estimateD(list( c(nSite, freqOrdr) ),
                               datatype = 'incidence_freq',
                               base = 'size', level = quotaN # , conf=NULL # 95% CI or none
                               )
    crRich <- crFull[crFull$order == 0, c('qD','qD.LCL','qD.UCL')]
    names(crRich) <- c('CRdiv','CRlow95','CRupr95')
    out <- cbind(out, crRich)
  }
  return(out)
}

# TODO return range sizes for all species in all subsamples
# taxDists <- function(dat, taxVar, idVar, sampIds,
#                      xy = NULL, omitSingles = FALSE){
#   if (omitSingles == TRUE){
#     freqs <- table( dat[,taxVar] )
#     singles <- names(freqs[freqs == 1])
#     rows2toss <- dat[,taxVar] %in% singles
#     dat <- dat[ !rows2toss, ]
#   }
#   # apply rangeSizer to all taxa
# }
