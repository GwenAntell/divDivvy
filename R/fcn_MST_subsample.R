# function to find all nested nearest-neighbour spatial subsamples
# modified from Roger Close's findAllNestedSpatialSubsamples.R
library(fields)
library(purrr) # for flatten(), but could write out longhand to reduce dependencies

# TODO discuss -
# * What should the output be? what format - list instead of tibble?
#   Return only a list of points within each cluster? ('cut the cake')
#   (Name each list element for seed if so)
# * Should branch cutoff be before or after max threshold reached?
# * Save only max cluster, not *everything* from min pts until max branch reached?
#   (e.g. list[i] = {1,2,3}, {1,2,3,4}, ... {1,2,3,...,n})

# for testing purposes
  # setwd("D:/Dropbox/OxfordPDRA")
  # datFls <- list.files('./data') %>% 
  #   grep('csv', ., value = TRUE)
  # pbdb_data <- paste0('data/', datFls[3]) %>%
  #   read.csv(header = TRUE, skip = 19)
  # lng_var <- "paleolng"
  # lat_var <- "paleolat"
  # min_n <- 3
  # max_branch <- 1000
  # i <- 1 # index ID of unique locations

# findAllNestedSpatialSubsamples # original function name
ptClustr <- function(pbdb_data, lng_var = 'paleolng', lat_var = 'paleolat',
                      # paleolngvar = "paleolng", paleolatvar = "paleolat", 
                     min_n = 3, max_branch = 1000 # max.subtree.size = 10000, 
                     # ID = "", # ID is bin name in RC script
                     # n_cores = n_cores # use this if adding option for parallel
                     ) {
  # RC: the bin name is fed to the function simply to add an ID string to spatial subsamples in output
  
	# set.seed(666)
  # TODO why set a seed? there aren't any random components
  # results identical between seed 6 and 666
  
  # ensure min number of occs, subset to unique locations, find GCD between all
	pcoords <- pbdb_data[,c(lng_var, lat_var)]
	colnames(pcoords) <- c('paleolng', 'paleolat')
	pbdb_data$pcoords <- paste(pbdb_data[[lng_var]], pbdb_data[[lat_var]], sep = "/")
	pcoords <- distinct(pcoords) %>% .[complete.cases(.), ] %>% as.data.frame()

	if (nrow(pcoords) < min_n) {
		stop('insufficient points for a cluster')
	  # return(NULL)
	}

	pcoords$ID <- paste("loc", 1:nrow(pcoords), sep = "")
	            # paste(ID, "loc", 1:nrow(pcoords), sep = "")
	pcoords <- pcoords[,c("ID","paleolng","paleolat")]
	nloc <- nrow(pcoords)

	gcd_matrix <- fields::rdist.earth(pcoords[,c("paleolng","paleolat")], miles = FALSE)
	# RC: One minor error in the function is that I didn't explicitly specify 
	# great circle distances computed by fields::rdist.earth() should be in km, 
	# so the sizes of spatial samples in the Science paper should be in miles
	colnames(gcd_matrix) <- rownames(gcd_matrix) <- pcoords$ID
	diag(gcd_matrix) <- NA

# find closest point, add it to set, and repeat until reach max tree branch
	ssFun <- function(i){  
		ss <- list()
		set <- rownames(gcd_matrix)[i]
		for (j in 1:(nloc)) {
	  #	find distance from each point in the set to all remaining points
			posscon <- gcd_matrix[which(pcoords$ID %in% set),
			                      -which(pcoords$ID %in% set), drop = F]
		# then retrieve name of closest next point
			nearest <- names(which.min(apply(posscon, MARGIN = 2, min)))
			set <- unique(c(set, nearest))
			
		#	set <- set[!is.na(set)] 
			if (set %>% is.na %>% any) stop('NA name in points set!')
			# flag an error if ever an NA value.
			# seems like there should never be, 
			# but RC must've added this check for a reason
			
			ss[[j]] <- set
			lngBranch <- max(gcd_matrix[which(pcoords$ID %in% set), 
			                            which(pcoords$ID %in% set)], na.rm = T)
			if (lngBranch > max_branch) break
		}
		return(ss)
	  } 
	ss <- lapply(1:nloc, ssFun)
	# if mc option specified with n_cores,
	# ss <- mclapply(seq_len(nrow(gcd_matrix)), ssFun,
	#	               mc.cores = n_cores, mc.preschedule = T, mc.cleanup = T, mc.silent = T)
	
	# flatten() removes 1 hierarchy level to squash list of lists into vector
	ss1 <- flatten(ss) %>% lapply(., sort) %>% unique
	
	# weed out any clusters with fewer than min specified points
	ss1 <- ss1[sapply(ss1, length) >= min_n]
	if (length(ss1) == 0){
	  # return(NULL)
	  stop('no set of points close enough to form a cluster')
	} 
	
	# back-calculate a df of lat-long coords for each cluster of names
	ss2 <- lapply(ss1, function(x) pcoords[which(pcoords$ID %in% x), ])
	names(ss2) <- paste("SS", 1:length(ss2), sep = "")
	           # paste(ID, "SS", 1:length(ss2), sep = "")

	# retrieve IDs of collections located at each cluster's pts
	# (note that multiple, unique collections can be at same point)
	coll_nos <- lapply(ss2, function(x) {
		temp <- paste(x$paleolng, x$paleolat, sep = "/")
		temp2 <- pbdb_data[which(pbdb_data$pcoords %in% temp), ]
		coll_nos <- unique(temp2$collection_no)
		return(coll_nos)
	})

	gcd_matrices <- lapply(ss2, function(x) {
		gcd_matrix[x$ID, x$ID]
	})
  
	# return everything that might be useful to query later
	ss_df <- tibble(SSID = names(ss2), SS = ss2, 
	                gcd_matrix = gcd_matrices, collection_nos = coll_nos)

	return(ss_df)

}
