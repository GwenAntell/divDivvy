library(terra)
library(sf)

# Use cell number/site ID and centroid coords as input.
# (Strictly, these don't have to be raster cell coords - could be coords of
# a predefined/clean set of locality points.)
# dat should contain unique site/cell IDs and locations
# (otherwise add step to subset pool to unique IDs)

# TODO option to pass on prj info to relevant function guts

# function of 1 subsample (1 try on 1 starting cell within 1 bin)
# take a subsample of sites/cells, without replacement, within buffered radius
cookie <- function(dat, siteId, xy, r, nSite, # prj,
                   seed, weight = FALSE # wasn't already an argument but should be
                 ){
  seedRow <- which(dat[, siteId] == seed)[1]
  # make sure coords are a 2-col, not-dataframe-class object to give to sf
  seedxy <- as.matrix( dat[seedRow, xy] )
  seedpt <- st_point(seedxy)
  # make sure to hard code the CRS for the buffer; sf can't infer lat-long
  seedsfc <- st_sfc(seedpt, crs = 'epsg:4326')
  buf <- st_buffer(seedsfc, dist = r*1000) 
  # also works:
  # rangemap::geobuffer_points(dat[1:5,xy], r*1000, wrap_antimeridian = TRUE)
  
  # split poly into multipolygon around antimeridian
  bufWrap <- st_wrap_dateline(buf, options = c("WRAPDATELINE=YES"))
  # examples apply st_union on this output, but seems unnecessary in tests here
  
  # find sites within radius of seed site/cell
  datSf <- st_as_sf(dat, coords = c('long','lat'),
                    crs = 'epsg:4326')
  poolBool <-  st_intersects(datSf, st_union(bufWrap), 
                             sparse = FALSE)
  # st_intersection would return geometries instead of logicals
  if (sum(poolBool) < nSite){
    stop('not enough pool cells to sample')
  }
  pool <- dat[poolBool, siteId]
  
  # assign sampling weights return ID names of rarefied sites/cells
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
    # it's ok weights don't sum to 1; sample() doesn't require this
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
