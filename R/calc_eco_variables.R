
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
#' Summarise the spatial extent and position of occurrence data, and
#' optionally estimate diversity and evenness
#'
#' \code{sdsumry} compiles a metadata about a sample or list of samples,
#' before or after spatial subsampling. The function counts the number
#' of unique spatial sites, collections (if requested), and taxa, and
#' calculates the centroid coordinates, latitudinal range (degrees),
#' great circle distance (km), and summed minimum spanning tree length (km)
#' for occurrences.
#'
#' If \code{quotaQ} is supplied, \code{sdsumry} rarefies richness at the
#' given coverage value and returns the point estimate of richness (Hill number 0)
#' and its 95% confidence interval as well as estimates of evenness (Pielou's J)
#' and sample coverage (Good's *u*). If \code{quotaN} is supplied,
#' \code{sdsumry} rarefies richness to the given number of occurrence counts
#' and returns the point estimate of richness and its 95% confidence interval.
#'
#' Coverage-based and classical rarefaction are both calculated with
#' \code{iNEXT::estimateD} internally. For details, such as how diversity
#' is extrapolated if a sample is insufficient to achieve a specified
#' rarefaction level, consult Chao and Jost (2012) and Hsieh *et al.* (2016).
#'
#' @param dat A \code{data.frame} or \code{matrix} containing taxon names,
#' coordinates, and any associated variables; or a list of such structures.
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
#' is supplied, remove the most common taxon prior to rarefaction. The `nTax`
#' and `evenness` returned are unaffected.
#'
#' @return A numeric vector of spatial and optional diversity metrics
#'
#' @export
#'
#' @references
#'
#' \insertRef{Chao2012}{divvy}
#'
#' \insertRef{Hsieh2016}{divvy}

sdsumry <- function(dat, taxVar, xy,
                    collections = NULL,
                    quotaQ = NULL, quotaN = NULL,
                    omitDom = FALSE){
  dfsumry <- function(df){
    out <- c()
    # n. collections (PBDB data) may be useful to know if rarefying
    if ( !is.null(collections) ){
      collSamp <- unique( df[,collections] )
      out <- cbind(out, 'nColl' = length(collSamp))
    }
    # comb out any duplicate occurrences of a taxon w/in single site.
    # do this after counting collections in case a single taxon
    # at a single site is recorded in multiple collections
    df <- unique( df[,c(taxVar, xy)] )

    df[,'site'] <- paste(df[, xy[1] ], df[, xy[2] ], sep = '/')
    siteIds <- unique( df[,'site'] )
    nSite <- length(siteIds)

    # run range size fcn as if all occs were from a single taxon
    # duplicate localities affect only occ count, nothing else
    sites <- unique( df[,xy] )
    spatSumry <- rangeSizer(coords = sites)
    out <- cbind(out, 'nSite'= nSite, spatSumry)

    # SCOR summed common occurrence rate (Hannisdal 2012)
    freqs <- table( df[,taxVar] )
    freqOrdr <- sort( as.numeric(freqs), decreasing = TRUE)
    # rate probability of detection under Poisson
    lambda <- -log(1 - freqOrdr/nSite) # natural log
    scor <- sum(lambda)

    # diversity metrics; optional coverage and/or classical rarefaction
    taxSamp <-  unique( df[,taxVar] )
    s <- length(taxSamp)
    out <- cbind(out, 'SCOR' = scor, 'nTax' = s) # raw count of taxa

    # Pielou's J evenness metric (calculate before omitting dom)
    pTax <- freqOrdr / sum(freqOrdr)
    h <- -sum(pTax * log(pTax))
    # vegan::diversity('shannon', base = exp(1))
    j <- h / log(s)

    if (omitDom == TRUE){
      freqOrdr <- freqOrdr[-1]
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
    out
  }
  if (class(dat) == 'list'){
    finL <- lapply(dat, dfsumry)
    fin <- data.frame(do.call(rbind, finL))
    # if original list has named elements e.g. latitudinal bands - keep them
    if ( ! is.null(names(dat)) ){
      fin$id <- names(dat)
    }
  }
  if (class(dat) %in% c('matrix', 'data.frame')){
    fin <- dfsumry(dat)
  }
  return(fin)
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
