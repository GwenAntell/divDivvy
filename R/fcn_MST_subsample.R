# function to find all nested nearest-neighbour spatial subsamples
# modified from Roger Close's findAllNestedSpatialSubsamples.R
library(sf)
# TODO build a find-seed fcn so all iters will be successful

# Changes to original (published) code -
# * Branch cutoff made *before* max threshold reached, not 1 point after
# * Save only max cluster, not *everything* from min pts until max branch reached
#   (e.g. list[i] = {1,2,3}, {1,2,3,4}, ... {1,2,3,...,n})
# * Start tree at random point and repeat for specified n iterations
#   (vs. exhaustively building at every possible start point)

# for testing purposes
  # datFls <- list.files('./data') %>%
  #   grep('csv', ., value = TRUE)
  # dat <- paste0('data/', datFls[1]) %>% read.csv
  # keepRows <- complete.cases(dat[,c('paleolng','paleolat','genus')])
  # dat <- dat[keepRows,]
  # xy <- c('paleolng', 'paleolat')
  # nMin <- 3
  # nSite <- 8 # if nSite supplied, nMin overridden
  # branchMax <- 1500 # in km
  # iter <- 500 # should be > 1 for sample() to work correctly
  # output <- 'locs'

ptClustr <- function(dat, xy, branchMax, nSite = NULL, iter, 
                     nMin = 3, output = 'locs') {
  # subset to unique locations and find dist matrix between all
	coords <- dat[,xy]
	dupes <- duplicated(coords)
	coords <- coords[ !dupes, ]
	coords <- coords[complete.cases(coords), ] 
	coords <- as.data.frame(coords) # in case data is given as a matrix
	# TODO only assign id if not already given? (e.g. cell name)
	coords$ID <- paste('loc', 1:nrow(coords), sep = '') 
	nLoc <- nrow(coords)
	if ( !is.null(nSite) ){
	  if (nLoc < nSite) stop('insufficient points for a subsample')
	} else {
	  if (nLoc < nMin) stop('insufficient points for a cluster')
	}
  
	coordSf <- st_as_sf(coords, coords = xy, crs = 'epsg:4326')
	gcdists <- st_distance(coordSf) # spherical distances (m) by default
	colnames(gcdists) <- rownames(gcdists) <- coords$ID
	diag(gcdists) <- NA
  
  # find closest point, add it to set, and repeat until reach max tree branch
	groupr <- function(id){
		clust <- id
		for (j in 1:(nLoc)) {
	  #	find distance from each point in the set to all remaining points
			posscon <- gcdists[ which(coords$ID %in% clust),
			                   -which(coords$ID %in% clust), drop = F]
		# then retrieve name of closest next point
			nearest <- names(which.min(apply(posscon, MARGIN = 2, min)))
			canddt <- unique(c(clust, nearest))
			# check if adding next pt puts cluster over size limit - if so don't add
			lngBranch <- max(gcdists[which(coords$ID %in% canddt), 
			                         which(coords$ID %in% canddt)], na.rm = T)
			branchMax <- set_units(branchMax, 'km')
			if (lngBranch > branchMax){
			  break
			}
			clust <- canddt
		}
		# TODO subsampling needs to happen BEFORE getting full df output
		# (which will contain duplicate coordinates from multitaxon sites)
		x <- xy[1]
		y <- xy[2]
		sampRows <- match(clust, coords$ID) # row location of sample pts in coord data
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
	}
	
	seeds <- sample(sample(coords$ID), iter, replace = TRUE)
	ss <- sapply(seeds, groupr, simplify = FALSE)
#	ss1 <- unique(ss) # duplicate subsample entities allowed, matches circ fcn
  
	# subsample trees to specified n sites
	if ( !is.null(nSite) ){
	  ss1 <- ss[sapply(ss, nrow) >= nSite]
	  ss1 <- lapply(ss1, function(clust){
	    clustRows <- sample(sample(1:nrow(clust)), nSite, replace = FALSE)
	    clust[clustRows,]
	  })
	} else {
	  ss1 <- ss[sapply(ss, nrow) >= nMin]
	}
	if (length(ss1) == 0){
	  stop('not enough close sites for any sample')
	} 
	
	return(ss1)
}
