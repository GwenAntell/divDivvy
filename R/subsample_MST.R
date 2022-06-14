# functions to find all nested nearest-neighbour spatial subsamples
# modified from Roger Close's findAllNestedSpatialSubsamples.R

# find closest point, add it to set, and repeat until reach max tree branch
# distMtrx arg = pairwise distances with units, rownames and colnames are pt IDs
groupr <- function(seed, sfPts, distMtrx, distMax){
  diag(distMtrx) <- NA
  distMax <- units::set_units(distMax, 'km')
  clust <- seed
  ptIds <- colnames(distMtrx)
  for (j in seq_along(ptIds)) { # alternative - while loop
    #	find distance from each point in the set to all remaining points
    distsOut <- distMtrx[ which(ptIds %in% clust),
                         -which(ptIds %in% clust), drop = F]
    # then retrieve name of closest next point
    mins <- apply(distsOut, MARGIN = 2, min)
    minMin <- min(mins)
    nearest <- names( which(mins == minMin) )
    # NB which.min() finds FIRST occurrence of min value, so instead:
    # hard-code to randomly pick in rare case of ties for min
    if (length(nearest) > 1){
      nearest <- sample(nearest, 1)
    }
    canddt <- c(clust, nearest)
    # check if adding next pt puts cluster over size limit - if so don't add
    diam <- max(distMtrx[which(ptIds %in% canddt),
                         which(ptIds %in% canddt)], na.rm = T)
    if (diam > distMax){
      break
    }
    clust <- canddt
  }
  return(clust)
}

#' Cluster localities within regions of nearest neighbours
#'
#' Spatially subsample a dataset based on minimum spanning trees connecting
#' points within regions of set extent, with optional rarefaction to a site quota.
#'
#' Lagomarcino and Miller (2012) developed an iterative approach of aggregating
#' localities to build clusters based on convex hulls, inspired by species-area
#' curve analysis (Scheiner 2003). Close et al. (2017, 2020) refined the approach and
#' changed the proximity metric from minimum convex hull area to minimum spanning
#' tree length. The present implementation adapts code from Close et al. (2020)
#' to add an option for site rarefaction after cluster construction and to grow
#' trees at random starting points \code{iter} number of times (instead of a
#' deterministic, exhaustive iteration at every unique location).
#'
#' The function takes a single location as a starting (seed) point; the seed
#' and its nearest neighbour initiate a spatial cluster. The distance between
#' the two points is the first branch in a minimum spanning tree for the cluster.
#' The location that has the shortest distance to any points already within the
#' cluster is grouped in next, and its distance (branch) is added to the sum
#' tree length. This iterative process continues until the largest distance
#' between any two points in the cluster would exceed \code{distMax} km.
#' In the very rare case multiple candidate points are tied for minimum distance
#' from the cluster, one point is selected at random as the next to include.
#' Any tree with fewer than \code{nMin} points is prohibited.
#'
#' In the case that \code{nSite} is supplied, \code{nMin} argument is ignored,
#' and any tree with fewer than \code{nSite} points is prohibited.
#' After building a tree as described above, a random set of \code{nSite} points
#' within the cluster is taken (without replacement). The returned output is
#' the coordinates of subsampled points if \code{output = 'locs'} or the subset
#' of \code{dat} associated with those coordinates if \code{output = 'full'}.
#' The \code{nSite} argument makes \code{clustr} comparable with \code{cookies}
#' in that it spatially standardises both extent and area/locality number.
#'
#' The performance of \code{clustr} is designed on the assumption \code{iter}
#' is much larger than the number of unique localities. Internal code first
#' calculates the full minimum spanning tree at every viable starting point
#' before it then samples those trees (i.e. resamples and optionally subsamples)
#' for the specified number of iterations. This sequence means the total
#' run-time increases only marginally even as \code{iter} increases greatly.
#' However, if there are a large number of sites, particularly a large number
#' of densely-spaced sites, the calculations will be slow even for a
#' small number of iterations.
#'
#' @inheritParams bandit
#' @param iter The number of spatial subsamples to return
#' @param distMax Numeric value for maximum diameter (km) allowed across
#' locations in a subsample
#' @param nMin Numeric value for the minimum number of sites to be included in
#' every returned subsample. If \code{nSite} supplied, \code{nMin} ignored.
#'
#' @seealso [cookies()]
#' @export
#' @importFrom Rdpack reprompt
#' @references
#'
#' \insertRef{Antell2020}{divvy}
#'
#' \insertRef{Close2017}{divvy}
#'
#' \insertRef{Close2020}{divvy}
#'
#' \insertRef{Lagomarcino2012}{divvy}
#'
#' \insertRef{Scheiner2003}{divvy}
#'

clustr <- function(dat, xy, iter, nSite = NULL, distMax, nMin = 3,
                   crs = 'epsg:4326', output = 'locs'){
  x <- xy[1]
  y <- xy[2]
	coords <- uniqify(dat[,xy], xy = xy)
	coords <- as.data.frame(coords) # in case data is given as a matrix
	nLoc <- nrow(coords)
	if ( !is.null(nSite) ){
	  if (nLoc < nSite) stop('insufficient points for a subsample')
	} else {
	  if (nLoc < nMin) stop('insufficient points for a cluster')
	}
	# great circle spherical distances for lon-lat coordinates (geodetic)
	# Euclidian distances for Cartesian coordinates
	coordSf <- sf::st_as_sf(coords, coords = xy, crs = crs)
	gcdists <- sf::st_distance(coordSf)
	coords$ID <- paste0('loc', 1:nrow(coords))
	colnames(gcdists) <- rownames(gcdists) <- coords$ID

	# build all possible trees, but return NULL if fewer pts than allowed
	posTrees <- sapply(coords$ID, function(seed){
	  tr <- groupr(seed, coordSf, gcdists, distMax)
	  if ( !is.null(nSite) ){
	    if (length(tr) >= nSite) tr
	  } else {
	    if (length(tr) >= nMin)  tr
	  }
	})
	availTrees <- Filter(Negate(is.null), posTrees)
	if (length(availTrees) == 0){
	  stop('not enough close sites for any sample')
	}

	smplTree <- function(){
	  # select one seed cell at random
	  seeds <- names(availTrees)
	  if (length(seeds) > 1){
	    seed <- sample(sample(seeds), 1)
	  } else {
	    # sample() fcn makes mistake if only 1 item to pick
	    seed <- seeds
	  }
	  # subsample tree from that seed
	  seedTr <- availTrees[seed][[1]]
	  if ( !is.null(nSite) ){
	    seedTr <- sample(sample(seedTr), nSite, replace = FALSE)
	  }
	  sampRows <- match(seedTr, coords$ID) # row location of sample pts in coord data
	  sampPtStrg <- paste(coords[sampRows, x], coords[sampRows, y], sep = '/')
	  datPtStrg  <- paste(dat[,x], dat[,y], sep = '/')
	  inSamp <- match(datPtStrg, sampPtStrg)
	  if (output == 'full'){
	    out <- dat[!is.na(inSamp), ]
	  } else {
	    if (output == 'locs'){
	      out <- coords[sampRows, xy]
	    } else {
	      stop('output argument must be one of c(\'full\', \'locs\')')
	    }
	  }
	  return(out)
	} # end function to take 1 subsample
	replicate(iter, smplTree(), simplify = FALSE)
}
