library(terra)
library(sf)
# TODO option to pass on prj info to relevant function guts

# return vector of cells that lie within buffer radius of given seed
# internal fcn for findPool and cookie
findPool <- function(seed, dat, siteId, xy, r, nSite # , prj
                     ){
  seedRow <- which(dat[, siteId] == seed)[1]
  # make sure coords are a 2-col, not-dataframe-class object to give to sf
  seedxy <- as.matrix( dat[seedRow, xy] )
  seedpt <- st_point(seedxy)
  # make sure to hard code the CRS for the buffer; sf can't infer lat-long
  seedsfc <- st_sfc(seedpt, crs = 'epsg:4326')
  buf <- st_buffer(seedsfc, dist = r*1000) 
  # split poly into multipolygon around antimeridian (patches 2020 bug)
  bufWrap <- st_wrap_dateline(buf, options = c("WRAPDATELINE=YES"))
  # identical results if st_union() applied to bufWrap
  # alternative way to construct buffer, but adding new pkg dependency:
  # rangemap::geobuffer_points(seedxy, r*1000, wrap_antimeridian = TRUE)
  
  # find sites within radius of seed site/cell
  datSf <- st_as_sf(dat, coords = c('long','lat'),
                    crs = 'epsg:4326')
  poolBool <-  st_intersects(datSf, bufWrap, sparse = FALSE)
  pool <- dat[poolBool, siteId]
  return(pool)
}

# function to try all possible starting pts (i.e. all occupied cells)
# save the ID of any cells that contain given pool size within buffer
findSeeds <- function(dat, siteId, xy, r, nSite, # prj
                      ){
  # count unique sites (not taxon occurences) relative to subsample quota
  dupes <- duplicated(dat[,siteId])
  dat <- dat[ !dupes, ]
  
  # test whether each occupied site/cell is viable for subsampling
  posSeeds <- dat[,siteId]
  isOK <- sapply(posSeeds, function(s){
    sPool <- findPool(s, dat, siteId, xy, r, nSite)
    n <- length(sPool)
    ifelse(n > nSite, s, NA)
  })
  # return IDs of only those seeds that are viable (pool > nSite)
  viable <- isOK[ !is.na(isOK)]
  return(viable)
}

# Use cell number/site ID and centroid coords as input.
# (Strictly, these don't have to be raster cell coords - could be coords of
# a predefined/clean set of locality points.)
# dat should contain unique site/cell IDs and locations
# (otherwise add step to subset pool to unique IDs)

# function of 1 subsample (1 try on 1 starting cell within 1 bin)
# take a subsample of sites/cells, without replacement, within buffered radius
cookies <- function(dat, siteId, xy, r, nSite, # prj,
                    seed, weight = FALSE){
  pool <- findPool(seed, dat, siteId, xy, r, nSite)
  if (length(pool) < nSite){
    stop('not enough pool cells to sample')
  }
  
  # assign sampling weights and return IDs of rarefied sites/cells
  if (weight){  
    samplIds <- sample(sample(pool), nSite, replace = FALSE) 
    # extra random sampling
    
  } else {
    # assign sampling weights inverse-square proportional to distance
    
    # remove seed from probabilistic sampling - include it manually
    pool <- pool[ !pool == seed]
    poolPts <- datSf[poolBool,]
    poolPts <- poolPts[ !poolPts$face == seed,]
    
    # squared inverse weight because inverse alone is too weak an effect
    gcdists <- st_distance(poolPts, seedsfc) # spherical distances by default
    wts <- sapply(gcdists, function(x) x^(-2))
    # sample() doesn't require wts sum = 1; identical results without rescaling
    samplIds <- c(seed, 
                  sample(sample(pool), nSite-1, prob = wts, replace = FALSE)
                  )
  }
  return(samplIds)
  
  # plot with global map to check it works correctly
    # library(maptools) # for plotting a simple world map
    # data("wrld_simpl")
    # plot(wrld_simpl, border="grey50")
    # plot(bufWrap, add=TRUE, col="red") # border = 'red', lwd=2
    # points(dat[,xy], pch=19, col="blue")
    # points(seedpt, pch=19, col="purple")
  # buffer area looks like inverse of what it should, but -
  # buffer does contain seed cell, so error must be 
  # in plot representation, not polygon geometry
    # st_intersects(seedsfc, bufWrap, sparse = FALSE) 
}
