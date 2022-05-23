#' Rarefy localities within latitudinal bands
#'
#' \code{bandit} subsamples spatial point data to a specified number of sites
#' within bins of equal latitude
#'
#' `bandit` rarefies the number of spatial sites within latitudinal ranges
#' of specified width. (Compare with `cookies` and `clustr`, which spatially
#' subsample to a specified total extent but without regard to position.)
#' Cases where it may be necessary to control for sampling site area by latitude
#' include characterisations of latitudinal diversity gradients (e.g. Marcot 2016)
#' or comparisons of environmental parameters that covary strongly with
#' latitude (e.g. diversity in reefal vs. non-reefal habitats). If there are
#' insufficient occurrences to draw a subsample in a given interval, that
#' interval is omitted from output. The latitudinal bounds of each subsample
#' are given by the name in the corresponding list element. Note that
#' the total surface area of the Earth within equal latitudinal increments
#' decreases from the equator towards the poles; `bandit` standardises only
#' the amount of sites/area encompassed by each subsample, not the total area
#' that could have been available for species to inhabit.
#'
#' @param dat A \code{data.frame} or \code{matrix} containing the
#' coordinate columns \code{xy} and any associated variables, e.g. taxon names.
#' @param xy A vector of two elements, specifying the name or numeric position
#' of the columns containing longitude and latitude coordinates, respectively.
#' @param width A numeric value for the resolution of latitude bins, in degrees.
#' @param iter The number of times to subsample localities within \strong{each}
#' latitudinal band.
#' @param nSite The quota of unique locations to include in each subsample.
#' @param centr Logical: should a bin center on and cover the equator
#' (\code{TRUE}) or should the equator mark the boundary between the
#' lowest-latitude northern and southern bins (\code{FALSE}, default)?
#' @param absLat Logical: should only the absolute values of latitude be
#' evaluated? If \code{absLat = TRUE}, \code{centr} argument is ignored.
#' @param output Whether the returned data should be a two-column matrix of
#' subsample site coordinates (\code{output = 'locs'}) or the subset of rows
#' from \code{dat} associated with those coordinates (\code{output = 'full'}).

#' @export
#'
#' @examples
#' data(bivalves)

#' # rasterise data into equal-area grid cells
#' library(icosa)
#' hgrid <- hexagrid( c(8,4) )
#' coords <- c('cellLng','cellLat')
#' faceIds <- locate(hgrid, bivalves[, c('paleolng','paleolat')] )
#' bivalves[, coords] <- pos(hgrid, faceIds)
#'
#' n <- 20
#' reps <- 100
#' set.seed(11)
#' # subsample within 10-degree bands of absolute latitude
#' bandAbs <- bandit(dat = bivalves, xy = coords,
#' iter = reps, nSite = n, output = 'full',
#' width = 10, absLat = TRUE
#' )
#' # head(bandAbs[[1]])
#' names(bandAbs)[1] # degree interval (absolute value)
#' #> "[0,10)"
#' unique(names(bandAbs)) # all intervals containing data
#' #> [1] "[0,10)" "[10,20)" "[30,40)" "[40,50)"
#'
#' # central latitude band spans equator, width = 20 degrees (as in Allen 2020)
#' # a finer-grain way to divide 180 degrees evenly into an odd number is width = 4
#' bandCent <- bandit(dat = bivalves, xy = coords,
#' iter = reps, nSite = n, output = 'full',
#' width = 20, centr = TRUE, absLat = FALSE
#' )
#' # head(bandCent[[1]])
#' names(bandCent)[1]
#' #> "[-10,10)"
#' unique(names(bandCent))
#' #> "[-10,10)" "[10,30)" "[30,50)"
#'
#' @references
#'
#' \insertRef{Allen2020}{divvy}
#'
#' \insertRef{Marcot2016}{divvy}

# TODO an option for equal-area latitudinal bands?
# TODO option for user-specified breaks

bandit <- function(dat, xy, iter, nSite, width,
                   centr = FALSE, absLat = FALSE,
                   output = 'locs'){
  x <- xy[1]
  y <- xy[2]
  coords <- uniqify(dat[,xy], xy = xy)
  lat <- coords[, y]
  if (absLat){
    lat <- abs(lat)
    brk <- seq(0, 90, by = width)
  } else {
    if (centr){
      # lowest-latitude band should straddle equator symmetrically
      brk <- c(seq(-width/2, -90, by = -width),
               seq( width/2,  90, by =  width))
      brk <- sort(brk)
    } else {
      brk <- seq(-90, 90, by = width)
    }
  }
  if (max(brk) != 90){
    warning(paste0('180 degrees latitude not evenly divisible by given width;',
                   ' northmost bound set below 90N'))
  }

  # end case of bin edge aligned at equator
  coords[,'band'] <- cut(lat, brk, right = FALSE) # labels arg
  # right = F needed so points at lat = 0 included, if lat is absolute val
  # (any point at N or S pole is problematic instead)
  bnds <- levels(coords[,'band'])

  # pick out bands with sufficient point density for subsampling
  bTally <- sapply(bnds, function(b) sum(coords[,'band'] %in% b) )
  bnds <- bnds[ bTally >= nSite ]
  if (length(bnds) < 1){
    stop('not enough close sites for any sample')
  }
  seeds <- sapply(bnds, replicate, n = iter)

  # rarefy sites within a band
  sampBnd <- function(b){
    bBool <- coords[,'band'] %in% b
    bDat <- coords[bBool,]
    sampRows <- sample(1:nrow(bDat), nSite, replace = FALSE)
    sampPtStrg <- paste(bDat[sampRows, x], bDat[sampRows, y], sep = '/')
    datPtStrg  <- paste(dat[,x], dat[,y], sep = '/')
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
