# function to find all nested nearest-neighbour spatial subsamples
# modified from Roger Close's findAllNestedSpatialSubsamples.R
library(sf)

# TODO discuss -
# * Should branch cutoff be before or after max threshold reached?
# * Save only max cluster, not *everything* from min pts until max branch reached?
#   (e.g. list[i] = {1,2,3}, {1,2,3,4}, ... {1,2,3,...,n})

# for testing purposes
  # datFls <- list.files('./data') %>%
  #   grep('csv', ., value = TRUE)
  # dat <- paste0('data/', datFls[3]) %>%
  #   read.csv(header = TRUE, skip = 19)
  # xy <- c('paleolng', 'paleolat')
  # nMin <- 3
  # branchMax <- 1500

ptClustr <- function(dat, xy = c('paleolng', 'paleolat'),
                     nMin = 3, branchMax = 1000 # max.subtree.size = 10000, 
                     ) {
  # subset to unique locations and find dist matrix between all
	coords <- dat[,xy]
	dupes <- duplicated(coords)
	coords <- coords[ !dupes, ]
	coords <- uniq[complete.cases(coords), ] 
	coords <- as.data.frame(coords) # in case data is given as a matrix
	# TODO only assign id if not already given? (e.g. cell name)
	coords$ID <- paste('loc', 1:nrow(coords), sep = "") 
	nLoc <- nrow(coords)
	if (nLoc < nMin) {
	  stop('insufficient points for a cluster')
	}
  
	coordSf <- st_as_sf(coords, coords = c(xVar, yVar),
	                    crs = 'epsg:4326')
	gcdists <- st_distance(coordSf) # spherical distances (m) by default
	colnames(gcdists) <- rownames(gcdists) <- coords$ID
	diag(gcdists) <- NA
  
  # find closest point, add it to set, and repeat until reach max tree branch
	groupr <- function(id){
#		ss <- list() # GSA edit: only return max cluster, not every clust from 3:max points
		clust <- id
		for (j in 1:(nLoc)) {
	  #	find distance from each point in the set to all remaining points
			posscon <- gcdists[ which(coords$ID %in% clust),
			                   -which(coords$ID %in% clust), drop = F]
		# then retrieve name of closest next point
			nearest <- names(which.min(apply(posscon, MARGIN = 2, min)))
			clust <- unique(c(clust, nearest))
			
#			ss[[j]] <- clust
			lngBranch <- max(gcdists[which(coords$ID %in% clust), 
			                         which(coords$ID %in% clust)], na.rm = T)
			units(branchMax) <- 'm'
			if (lngBranch > branchMax * 1000) break
		}
		sort(clust) # return(ss)
	}
	
	ss <- sapply(coords$ID, groupr)
	ss1 <- unique(ss)
  
	# weed out any clusters with fewer than min specified points
	ss1 <- ss1[sapply(ss1, length) >= nMin]
	if (length(ss1) == 0){
	  stop('not enough close sites for any sample')
	} 
	
	return(ss1)
}
