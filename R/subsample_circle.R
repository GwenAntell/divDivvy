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
  # don't use sapply or object will condense from list to matrix
  # in the special case all pool lengths are equal:
  posPools <- lapply(posSeeds, function(s){
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
#' the antemeridian (180 degrees longitude) are wrapped as a multipolygon
#' to prevent artificial truncation. After standardising radial extent, sites
#' are drawn within the circular extent until a quota of \code{nSite} is met.
#' Sites are sampled without replacement, so a location is used as a seed point
#' only if it is within \code{r} km distance of at least \code{nSite-1} locations.
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
#' @param r Numeric value for the radius (km) defining the circular extent
#' of each spatial subsample.
#' @param weight Whether sites within the subsample radius should be drawn
#' at random (\code{weight = FALSE}, default) or with probability inversely
#' proportional to the square of their distance from the centre of the
#' subsample region (\code{weight = TRUE}).
#'
#' @returns A list of length \code{iter}. Each list element is a
#' \code{data.frame} or \code{matrix} (matching the class of \code{dat})
#' with \code{nSite} observations. If \code{output = 'locs'}
#' (default), only the coordinates of subsampling locations are returned.
#' If \code{output = 'full'}, all \code{dat} columns are returned for the
#' rows associated with the subsampled locations.
#'
#' If \code{weight = TRUE}, the first observation in each returned subsample
#' \code{data.frame} corresponds to the seed point. If \code{weight = FALSE},
#' observations are listed in the random order of which they were drawn.
#'
#' @examples
#' # generate occurrences: 10 lat-long points in modern Australia
#' n <- 10
#' x <- seq(from = 140, to = 145, length.out = n)
#' y <- seq(from = -20, to = -25, length.out = n)
#' pts <- data.frame(x, y)
#'
#' # sample 5 sets of 3 occurrences within 200km radius
#' cookies(dat = pts, xy = 1:2, iter = 5,
#'         nSite = 3, r = 200)
#'
#' @seealso [clustr()]
#' @export
#' @references
#'
#' \insertRef{Antell2020}{divvy}
cookies <- function(dat, xy, iter, nSite, r, weight = FALSE,
                    crs = 'epsg:4326', output = 'locs'){
  coords <- uniqify(dat, xy) |> as.data.frame()
  coords$id <- paste0('loc', 1:nrow(coords))

  # this is the rate-limiting step (very slow), but overall
  # it's more efficient to construct all spatial buffers here at start
  # and not repeat the calculations anywhere later!
  allPools <- findSeeds(coords, 'id', xy, r, nSite, crs)
  if (length(allPools) < 1){
    stop('not enough close sites for any subsample')
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
      poolBool <- coords[, 'id'] %in% pool
      poolPts <- datSf[poolBool,]

      # squared inverse weight because inverse alone is too weak an effect
      # great circle spherical distances for lon-lat coordinates (geodetic)
      # Euclidian distances for Cartesian coordinates
      seedRow <- which(coords[, 'id'] == seed)[1]
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
    coordRows <- match(samplIds, coords$id) # row location of sample pts in coord data
    coordLocs <- coords[coordRows, xy]

    if (output == 'full'){
      x <- xy[1]
      y <- xy[2]
      sampPtStrg <- paste(coordLocs[, x], coordLocs[, y], sep = '/')
      datPtStrg  <- paste(      dat[, x],       dat[, y], sep = '/')
      inSamp <- match(datPtStrg, sampPtStrg)
      out <- dat[ !is.na(inSamp), ]
    } else {
      if (output == 'locs'){
        out <- coordLocs
      } else {
        stop('output argument must be one of c(\'full\', \'locs\')')
      }
    } # end output formatting
    return(out)
  }
  replicate(iter, cookie(), simplify = FALSE)
}
