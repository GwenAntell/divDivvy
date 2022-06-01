# return vector of cells that lie within buffer radius of given seed
findPool <- function(seed, dat, siteId, xy, r, nSite, crs = 'epsg:4326'
                     ){
  datSf <- sf::st_as_sf(dat, coords = xy, crs = crs)
  seedRow <- which(dat[, siteId] == seed)[1]
  seedpt <- datSf[seedRow, ]
  # buffer will be more accurate if projected,
  # but wrapping around antimeridian requires lat-long coordinates
  r <- units::set_units(r, 'km')
  buf <- sf::st_buffer(seedpt, dist = r)
  if (crs != 'epsg:4326'){
    buf   <- sf::st_transform(buf,   crs = 'epsg:4326')
    datSf <- sf::st_transform(datSf, crs = 'epsg:4326')
  }
  bufWrap <- sf::st_wrap_dateline(buf, options = c("WRAPDATELINE=YES"))

  # find sites within radius of seed site/cell
  poolBool <- sf::st_intersects(datSf, bufWrap, sparse = FALSE)
  pool <- dat[poolBool, siteId]
  return(pool)
}

# function to try all possible starting pts (i.e. all occupied cells)
# save the ID of any cells that contain given pool size within buffer
findSeeds <- function(dat, siteId, xy, r, nSite, crs = 'epsg:4326'
                      ){
  # test whether each occupied site/cell is viable for subsampling
  posSeeds <- dat[,siteId]
  posPools <- sapply(posSeeds, function(s){
    sPool <- findPool(s, dat, siteId, xy, r, nSite, crs)
    n <- length(sPool)
    if (n >= nSite)
      sPool
  })
  # return pool site/cell IDs for each viable seed point
  # same overall list structure as cookies outputs; names = seed IDs
  names(posPools) <- posSeeds
  Filter(Negate(is.null), posPools)
}


#' Rarefy localities within circular regions of standard area
#'
#' Spatially subsample a dataset to produce samples of standard area and extent.
#'
#' The function takes a single location as a starting (seed) point and
#' circumscribes a buffer of \code{r} km around it. Buffer circles that span
#' the antemeridian (180 deg longitude) are wrapped as a multipolygon
#' to prevent artificial truncation. After standardising radial extent, sites
#' are drawn within the circular extent until a quota of \code{nSite}.
#' Sites are sampled without replacement, so a location is used as a seed point
#' only if it is within \code{r} km distance of at least \code{nSite} locations.
#' The method is introduced in Antell et al. (2020) and described in
#' detail in Methods S1 therein.
#'
#' The probability of drawing each site within the standardised extent is
#' either equal (\code{weight = FALSE}) or proportional to the inverse-square
#' of its distance from the seed point (\code{weight = TRUE}), which clusters
#' subsample locations more tightly.
#'
#' For geodetic coordinates (latitude-longitude), distances are calculated along
#' great circle arcs. For Cartesian coordinates, distances are calculated in
#' Euclidian space, in units associated with the projection CRS (e.g. metres).
#'
#' @inheritParams clustr
#' @param siteId The name or numeric position of the column in \code{dat}
#' containing identifiers for unique spatial sites, e.g. raster cell names.
#' @param xy A vector of two elements, specifying the name or numeric position
#' of columns in \code{dat} containing coordinates, e.g. longitude and latitude.
#' Coordinates for the same site ID should be identical, and where IDs are
#' raster cells the coordinates are usually expected to be cell centroids.
#' @param r Numeric value for the radius (km) defining the circular extent
#' of each spatial subsample.
#' @param weight Whether sites within the subsample radius should be drawn
#' at random (\code{weight = FALSE}) or with probability inversely proportional
#' to the square of their distance from the centre of the subsample region.
#' @param crs Coordinate reference system as a GDAL text string, numeric EPSG
#' code, or object of class `crs`. Default is latitude-longitude (EPSG:4326).
#' @param output Whether the returned data should be a two-column matrix of
#' subsample site coordinates (\code{output = 'locs'}), with row names the
#' site IDs, or the subset of rows from \code{dat} associated with those
#' coordinates (\code{output = 'full'}).
#'
#' @seealso [clustr()]
#' @export
#' @references
#'
#' \insertRef{Antell2020}{divvy}
cookies <- function(dat, xy, iter, nSite, siteId, r, weight = FALSE,
                    crs = 'epsg:4326', output = 'locs'){
  locDat <- dat[, c(xy, siteId)]
  coords <- uniqify(locDat, siteId = siteId)

  # this is the rate-limiting step (v slow), but overall
  # it's most efficient to construct all spatial buffers here at start
  # and not repeat the calculations anywhere later!
  allPools <- findSeeds(coords, siteId, xy, r, nSite, crs)
  if (length(allPools) < 1){
    stop('not enough close sites for any sample')
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
      # convert to spatial features for distance calculations
      datSf <- sf::st_as_sf(coords, coords = xy, crs = crs)

      # remove seed from probabilistic sampling - include it manually
      # (otherwise inverse distance will divide by zero)
      pool <- pool[ !pool == seed]
      poolBool <- coords[,siteId] %in% pool
      poolPts <- datSf[poolBool,]

      # squared inverse weight because inverse alone is too weak an effect
      # great circle spherical distances for lon-lat coordinates (geodetic)
      # Euclidian distances for Cartesian coordinates
      seedRow <- which(coords[, siteId] == seed)[1]
      seedPt <- datSf[seedRow,]
      gcdists <- sf::st_distance(poolPts, seedPt)
      wts <- sapply(gcdists, function(x) x^(-2))
      # sample() doesn't require wts sum = 1; identical results without rescaling
      samplIds <- c(seed,
                    sample(sample(pool), nSite-1, prob = wts, replace = FALSE)
      )
    } else {
      samplIds <- sample(sample(pool), nSite, replace = FALSE)
    } # end site rarefaction
    if (output == 'full'){
      inSamp <- match(dat[, siteId], samplIds)
      out <- dat[!is.na(inSamp), ]
    } else {
      if (output == 'locs'){
        inSamp <- match(samplIds, coords[, siteId])
        out <- coords[inSamp, xy]
        rownames(out) <- samplIds
      } else {
        stop('output argument must be one of c(\'full\', \'locs\')')
      }
    } # end output formatting
    return(out)
  }
  replicate(iter, cookie(), simplify = FALSE)
}
