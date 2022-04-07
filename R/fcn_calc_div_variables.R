library(iNEXT)
library(units)

# for testing purposes
# dat <- sites
# taxVar <- 'genus'
# idVar <- 'face' # rename to 'siteVar' or 'sampVar' ?
# sampIds <- circs[[1]]
# collections <- 'collection_no' # remove as too problematic a variable?
# xy <- c('long','lat') # c('paleolng','paleolat') # if using original, vector coords
# quotaQ <- 0.2
# quotaN <- 10
# omitDom <- TRUE

# dat should be list of occs for a single taxon, which are presumed unique
rangeSizer <- function(dat, xy = NULL){
  out <- cbind('nOcc' = nrow(dat))
  if ( !is.null(xy) ){
    coords <- dat[,xy]
    latDiff <- max(coords[,2]) - min(coords[,2])
    latRng <- abs(latDiff)
    pts <- st_as_sf(coords, coords = xy, crs = 'epsg:4326')
    ptsGrp <- st_union(pts)
    cntr <- unlist( st_centroid(ptsGrp) )
    gcdists <- st_distance(pts) # returns units-class object (m)
    gcdists <- set_units(gcdists, 'km')
    gcMax <- max(gcdists) # / 1000 
    mst <- spantree(gcdists) # drop_units(gcdists)
    agg <- sum(mst$dist) # / 1000 
    out <- cbind(out,
                 'centroidLng' = cntr[1],
                 'centroidLat' = cntr[2],
                 'latRange' = latRng,
                 'greatCircDist' = gcMax,
                 'minSpanTree' = agg
    )
  }
  return(out)
}

# collections arg = name of col containing collection id, e.g. 'collection_no'
# coords arg = name of latlong cols to calculate lat range, MST, and midpoint
# quotaN = if given, perform classical rarefaction w specified quota (n occs)
# quotaQ = if given, perform SQS w specified coverage/quorum level
# omitX = removes single/dominant species from ALL metrics

# return: great circle distance and summed tree length is returned in km
# if too few taxon occs to achieve specified rarefaction level, div is extrap

sampMeta <- function(dat, taxVar, idVar, sampIds,
                     collections = NULL, xy = NULL,
                     quotaQ = NULL, quotaN = NULL, 
                     omitDom = FALSE){
  smpld <- dat[, idVar] %in% sampIds
  datSamp <- dat[smpld,]
  out <- c()
  # n. collections (PBDB data) may be useful to know if rarefying
  if ( !is.null(collections) ){
    collSamp <- unique( datSamp[,collections] )
    out <- cbind(out, 'nColl' = length(collSamp))
  }
  # comb out any duplicate occurrences of a taxon w/in single site.
  # do this after counting collections in case a single taxon
  # at a single site is recorded in multiple collections
  if (is.null(xy)){
    datSamp <- unique( datSamp[,c(taxVar, idVar)] )
  } else {
    datSamp <- unique( datSamp[,c(taxVar, idVar, xy)] )
  }
  
  siteSamp <- unique( datSamp[,idVar])
  nSite <- length(siteSamp)
  out <- cbind(out, 'nSite'= nSite)

  # run range size fcn as if all occs were from a single taxon
  # duplicate localities affect only occ count, nothing else
  spatMeta <- rangeSizer(datSamp, xy = xy)
  out <- cbind(out, spatMeta)
  
  # diversity metrics; optional coverage and/or classical rarefaction
  taxSamp <-  unique( datSamp[,taxVar] )
  out <- cbind(out, 'nTax' = length(taxSamp)) # raw count of taxa
  freqs <- table( datSamp[,taxVar] )
  freqOrdr <- sort( as.numeric(freqs), decreasing = TRUE)
  if (omitDom == TRUE){
    dom <- names( which.max(freqs) )
    domRows <- datSamp[,taxVar] == dom
    datSamp <- datSamp[ !domRows, ]
  }
  if ( !is.null(quotaQ) ){
    sqsFull <- estimateD(list( c(nSite, freqOrdr) ), 
                         datatype = 'incidence_freq',
                         base = 'coverage', level = quotaQ # , conf=NULL # 95% CI or none
                         )
    # q0 = richness, q1 = exp of Shannon's entropy index, 
    # q2 = inverse of Simpson's concentration index
    sqsRich <- sqsFull[sqsFull$order == 0, c('SC','qD','qD.LCL','qD.UCL')]
    names(sqsRich) <- c('u','SQSdiv','SQSlow95','SQSupr95')
    # return sample coverage, richness estimate and 95% CI
    out <- cbind(out, sqsRich)
  }
  if ( !is.null(quotaN) ){
    crFull <- estimateD(list( c(nSite, freqOrdr) ), 
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
taxDists <- function(dat, taxVar, idVar, sampIds,
                     xy = NULL, omitSingles = FALSE){
  if (omitSingles == TRUE){
    freqs <- table( dat[,taxVar] )
    singles <- names(freqs[freqs == 1])
    rows2toss <- dat[,taxVar] %in% singles
    dat <- dat[ !rows2toss, ]
  }
  # apply rangeSizer to all taxa
}
