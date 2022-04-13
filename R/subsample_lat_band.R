#' Rarefy localities within latitudinal bands
#'
#' \code{bandit} subsamples spatial point data to a specified number of sites
#' within bins of equal latitude
#'
#' @param dat A \code{data.frame} or \code{matrix} containing the
#' coordinate columns \code{xy} and any associated variables, e.g. taxon names.
#' @param xy A vector of two elements, specifying the name or numeric position
#' of the columns containing longitude and latitude coordinates.
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
#' from \code{dat} associated with those coordinates (\code{output = 'all'}).

#' @export

# TODO an option for equal-area latitudinal bands?

# output one of 'all' or 'locs'
bandit <- function(dat, xy, width, iter, nSite,
                   centr = FALSE, absLat = FALSE,
                   output = 'locs'){
  x <- xy[1]
  y <- xy[2]
  dupes <- duplicated(dat[,xy])
  coords <- dat[ !dupes, ]
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
  # if (max(brk) != 90){
  #   stop('latitude not evenly divisible by width')
  # }

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
