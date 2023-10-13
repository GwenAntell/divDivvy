#' Rarefy localities within latitudinal bands
#'
#' `bandit` subsamples spatial point data to a specified number of sites
#' within bins of equal latitude
#'
#' `bandit` rarefies the number of spatial sites within latitudinal ranges
#' of specified bin width. (Compare with `cookies` and `clustr`, which spatially
#' subsample to a specified extent without regard to latitudinal position.)
#' Cases where it could be appropriate to control for latitudinal spread of localities
#' include characterisations of latitudinal diversity gradients (e.g. Marcot 2016)
#' or comparisons of ecosystem parameters that covary strongly with
#' latitude (e.g. diversity in reefal vs. non-reefal habitats). Note that
#' the total surface area of the Earth within equal-latitudinal increments
#' decreases from the equator towards the poles; `bandit` standardises only
#' the amount of sites/area encompassed by each subsample, not the total area
#' that could have been available for species to inhabit.
#'
#' As with all `divvy` subsampling functions, sites within a given
#' regional/latitudinal subsample are selected without replacement.
#'
#' To calculate an integer number of degrees into which a given latitudinal
#' range divides evenly, the `palaeoverse` package (v 1.2.1) provides the
#' `lat_bins` function with argument `fit = TRUE`.
#'
#' @param dat A `data.frame` or `matrix` containing the coordinate
#' columns `xy` and any associated variables, e.g. taxon names.
#' @param xy A vector of two elements, specifying the name or numeric position
#' of columns in `dat` containing coordinates, e.g. longitude and latitude.
#' Coordinates for any shared sampling sites should be identical, and where sites
#' are raster cells, coordinates are usually expected to be cell centroids.
#' @param crs Coordinate reference system as a GDAL text string, EPSG code,
#' or object of class `crs`. Default is latitude-longitude (`EPSG:4326`).
#' @param bin A positive numeric value for latitudinal band width, in degrees.
#' @param iter The number of times to subsample localities within \strong{each}
#' latitudinal band.
#' @param nSite The quota of unique locations to include in each subsample.
#' @param centr Logical: should a bin center on and cover the equator
#' (`TRUE`) or should the equator mark the boundary between the
#' lowest-latitude northern and southern bins (`FALSE`, default)?
#' Ignored if `absLat = TRUE`.
#' @param absLat Logical: should only the absolute values of latitude be
#' evaluated? If `absLat = TRUE`, `centr` argument is ignored.
#' @param maxN Optional argument to specify the northmost limit for subsampling,
#' if less than 90 degrees.
#' @param maxS Optional argument to specify the southmost limit for subsampling,
#' if not -90 degrees. Should be a negative value if in the southern hemisphere.
#' @param output Whether the returned data should be two columns of
#' subsample site coordinates (`output = 'locs'`) or the subset of rows
#' from `dat` associated with those coordinates (`output = 'full'`).
#'
#' @return A list of subsamples, each a `data.frame` containing
#' coordinates of subsampled localities (if `output = 'locs'`)
#' or the subset of occurrences from `dat` associated with those coordinates
#' (if `output = 'full'`). The latitudinal bounds of each subsample
#' are specified by its name in the list. If there are too few localities
#' in a given interval to draw a subsample, that interval is omitted from output.
#'
#' @seealso [cookies()]
#' @seealso [clustr()]
#' @export
#'
#' @examples
#' # load bivalve occurrences to rasterise
#' library(terra)
#' data(bivalves)
#'
#' # initialise Equal Earth projected coordinates
#' rWorld <- rast()
#' prj <- 'EPSG:8857'
#' rPrj <- project(rWorld, prj, res = 200000) # 200,000m is approximately 2 degrees
#'
#' # coordinate column names for the current and target coordinate reference system
#' xyCartes <- c('paleolng','paleolat')
#' xyCell   <- c('centroidX','centroidY')
#'
#' # project occurrences and retrieve cell centroids in new coordinate system
#' llOccs <- vect(bivalves, geom = xyCartes, crs = 'epsg:4326')
#' prjOccs <- project(llOccs, prj)
#' cellIds <- cells(rPrj, prjOccs)[,'cell']
#' bivalves[, xyCell] <- xyFromCell(rPrj, cellIds)
#'
#' # subsample 20 equal-area sites within 10-degree bands of absolute latitude
#' n <- 20
#' reps <- 100
#' set.seed(11)
#' bandAbs <- bandit(dat = bivalves, xy = xyCell,
#'                   iter = reps, nSite = n, output = 'full',
#'                   bin = 10, absLat = TRUE,
#'                   crs = prj
#' )
#' head(bandAbs[[1]]) # inspect first subsample
#' names(bandAbs)[1] # degree interval (absolute value) of first subsample
#' #> [1] "[10,20)"
#' unique(names(bandAbs)) # all intervals containing sufficient data
#' #> [1] "[10,20)" "[20,30)" "[30,40)" "[40,50)"
#' # note insufficient coverage to subsample at equator or above 50 degrees
#'
#' # subsample 20-degree bands, where central band spans the equator
#' # (-10 S to 10 N latitude), as in Allen et al. (2020)
#' # (An alternative, finer-grain way to divide 180 degrees evenly into an
#' # odd number of bands would be to set 'bin' = 4.)
#' bandCent <- bandit(dat = bivalves, xy = xyCell,
#'                    iter = reps, nSite = n, output = 'full',
#'                    bin = 20, centr = TRUE, absLat = FALSE,
#'                    crs = prj
#' )
#' unique(names(bandCent)) # all intervals containing sufficient data
#' #> [1] "[-50,-30)" "[10,30)" "[30,50)"
#'
#' @references
#'
#' \insertRef{Allen2020}{divvy}
#'
#' \insertRef{Marcot2016}{divvy}

bandit <- function(dat, xy, iter, nSite, bin,
                   centr = FALSE, absLat = FALSE,
                   maxN = 90, maxS = -90,
                   crs = 'epsg:4326', output = 'locs'){
  coords <- uniqify(dat[,xy], xy = xy)
  coords <- as.data.frame(coords) # in case data is given as a matrix
  sfCoords <- sf::st_as_sf(coords, coords = xy, crs = crs)
  if ( ! sf::st_is_longlat(sfCoords) ){
    ll <- sf::sf_project(from = crs, to = 'epsg:4326', coords,
                         keep = TRUE, warn = TRUE)
    coords <- cbind(coords, 'long' = ll[,1], 'lat' = ll[,2])
    latCol <- 'lat'
  } else {
    latCol <- xy[2]
  }

  lat <- coords[, latCol]
  if (absLat){
    lat <- abs(lat)
    maxAbs <- c(maxN, maxS) |> abs() |> max()
    brk <- seq(0, maxAbs, by = bin)

    if (max(brk) != maxAbs){
      warning(paste0('latitudinal range not evenly divisible by given bin width:',
                     ' pole-most bound set less than given maximum'))
    }
  } else {
    if (centr){
      # lowest-latitude band should straddle equator symmetrically
      brk <- c(seq(-bin/2, maxS, by = -bin),
               seq( bin/2,  maxN, by =  bin))
      brk <- sort(brk)
    } else {
      brk <- seq(maxS, maxN, by = bin)
    }

    if (max(brk) != maxN){
      warning(paste0('latitudinal range not evenly divisible by given bin width:',
                     ' northmost bound set less than given maximum'))
    }
    if (min(brk) != maxS){
      warning(paste0('latitudinal range not evenly divisible by given bin width:',
                     ' southmost bound set less than given maximum'))
    }
  } # end case of signed lat values

  # end case of bin edge aligned at equator
  coords[,'band'] <- cut(lat, brk, right = FALSE) # labels arg
  # right = F needed so points at lat = 0 included, if lat is absolute val
  # (any point at N or S pole is problematic instead)
  bnds <- levels(coords[,'band'])

  # pick out bands with sufficient point density for subsampling
  bTally <- table(coords[,'band']) |> as.numeric()
  bnds <- bnds[ bTally >= nSite ]
  if (length(bnds) < 1){
    stop('not enough close sites for any subsample')
  }
  seeds <- sapply(bnds, replicate, n = iter)

  # rarefy sites within a band
  x <- xy[1]
  y <- xy[2]
  datPtStrg  <- paste(dat[,x], dat[,y], sep = '/')
  sampBnd <- function(b){
    bBool <- coords[,'band'] %in% b
    bDat <- coords[bBool,]
    sampRows <- sample(1:nrow(bDat), nSite, replace = FALSE)
    sampPtStrg <- paste(bDat[sampRows, x], bDat[sampRows, y], sep = '/')
    inSamp <- match(datPtStrg, sampPtStrg)

    if (output == 'full'){
      out <- dat[!is.na(inSamp), ]
    } else {
      if (output == 'locs'){
        out <- bDat[sampRows, xy]
      } else {
        stop('output argument must be one of c(\'full\', \'locs\')')
      }
    }
    return(out)
  }
  sapply(seeds, sampBnd, simplify = FALSE)
}
