# function to find all nested nearest-neighbour spatial subsamples
# modified from Roger Close's findAllNestedSpatialSubsamples.R

# Changes to original (published) code -
# * Branch cutoff made *before* max threshold reached, not 1 point after
# * Save only max cluster, not *everything* from min pts until max branch reached
#   (e.g. list[i] = {1,2,3}, {1,2,3,4}, ... {1,2,3,...,n})
# * Start tree at random point and repeat for specified n iterations
#   (vs. exhaustively and uniquely building at every point)

# for a lot of paleo datasets, there will be fewer unique sites
# than desired number of iterations, so more efficient to build every possible tree
# up front and then (re)sample them n iterations, instead of
# recalculating a new tree each iteration

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
    nearest <- names(which.min(apply(distsOut, MARGIN = 2, min)))
    canddt <- c(clust, nearest) # unique
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
#' curve analysis (Scheiner 2003). Close et al. (2020) refined the approach and
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
#' \insertRef{Antell2020}{divDivvy}
#'
#' \insertRef{Close2020}{divDivvy}
#'
#' \insertRef{Lagomarcino2012}{divDivvy}
#'
#' \insertRef{Scheiner2003}{divDivvy}
#'

clustr <- function(dat, xy, distMax, nSite = NULL, iter,
                   nMin = 3, output = 'locs') {
  # subset to unique locations and find dist matrix between all
  x <- xy[1]
  y <- xy[2]
  coords <- dat[,xy]
	dupes <- duplicated(coords)
	coords <- coords[ !dupes, ]
	coords <- coords[stats::complete.cases(coords), ]
	coords <- as.data.frame(coords) # in case data is given as a matrix
	nLoc <- nrow(coords)
	if ( !is.null(nSite) ){
	  if (nLoc < nSite) stop('insufficient points for a subsample')
	} else {
	  if (nLoc < nMin) stop('insufficient points for a cluster')
	}
	coordSf <- sf::st_as_sf(coords, coords = xy, crs = 'epsg:4326')
	gcdists <- sf::st_distance(coordSf) # spherical distances (m) by default
	coords$ID <- paste0('loc', seq_along(coords))
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
	  # subsampling needs to happen BEFORE getting full df output
	  # (which can contain duplicate coordinates from multitaxon sites)
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
