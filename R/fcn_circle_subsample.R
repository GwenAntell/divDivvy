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
  # identical results if st_union() applied to bufWrap.
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
findSeeds <- function(dat, siteId, xy, r, nSite # , prj
                      ){
  # count unique sites (not taxon occurences) relative to subsample quota
  dupes <- duplicated(dat[,siteId])
  dat <- dat[ !dupes, ]
  
  # test whether each occupied site/cell is viable for subsampling
  posSeeds <- dat[,siteId]
  posPools <- sapply(posSeeds, function(s){
    sPool <- findPool(s, dat, siteId, xy, r, nSite)
    n <- length(sPool)
    if (n > nSite)
      sPool
  })
  # return pool site/cell IDs for each viable seed point
  # same overall list structure as cookies outputs; names = seed IDs
  Filter(Negate(is.null), posPools)
  # posPools[!sapply(posPools, is.null)] # equivalent, base syntax
  
}

# Use cell number/site ID and centroid coords as input.
# (Strictly, these don't have to be raster cell coords - could be coords of
# a predefined/clean set of locality points.)
# dat should contain unique site/cell IDs and locations
# (otherwise add step to subset pool to unique IDs)

cookies <- function(dat, siteId, xy, r, nSite, # prj,
                    iter, weight = FALSE){
  # this is the rate-limiting step (v slow), but overall
  # it's most efficient to construct all spatial buffers here at start
  # and not repeat the calculations anywhere later!
  allPools <- findSeeds(dat, siteId, xy, r, nSite)
  if (length(allPools) < 1){
    stop('not enough close sites for any sample')
  }
  # convert to spatial features for distance calculations later
  if (weight){
    datSf <- st_as_sf(dat, coords = c('long','lat'),
                      crs = 'epsg:4326')
  }
  
  # takes a subsample of sites/cells, w/o replacement, w/in buffered radius
  cookie <- function(){
    # select one seed cell at random
    seeds <- names(allPools)
    if (length(seeds) > 1){
      seed <- sample(sample(seeds), 1) 
    } else { 
      # sample() fcn makes mistake if only 1 item to pick
      seed <- seeds
    }
    pool <- allPools[seed][[1]]
    
    if (weight){  
      # remove seed from probabilistic sampling - include it manually
      # (otherwise inverse distance will divide by zero)
      pool <- pool[ !pool == seed]
      poolBool <- dat[,siteId] %in% pool
      poolPts <- datSf[poolBool,] 
      
      # squared inverse weight because inverse alone is too weak an effect
      seedRow <- which(dat[, siteId] == seed)[1]
      seedPt <- datSf[seedRow,]
      gcdists <- st_distance(poolPts, seedPt) # spherical distances by default
      wts <- sapply(gcdists, function(x) x^(-2))
      # sample() doesn't require wts sum = 1; identical results without rescaling
      samplIds <- c(seed, 
                    sample(sample(pool), nSite-1, prob = wts, replace = FALSE)
      )
    } else {
      samplIds <- sample(sample(pool), nSite, replace = FALSE) 
    }
    return(samplIds) # IDs of rarefied sites/cells
  }
  replicate(iter, cookie(), simplify = FALSE)
}
