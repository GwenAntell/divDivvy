# datFls <- list.files('./data') %>%
#   grep('csv', ., value = TRUE)
# df <- paste0('data/', datFls[1]) %>% read.csv
# #  read.csv(header = TRUE, skip = 19) # if metadata included in download file
# keepRows <- complete.cases(df[,c('paleolng','paleolat','genus')])
# df <- df[keepRows,]
# 
# dat <- df[,c('genus','paleolng','paleolat')]
# xy <- c('paleolng','paleolat')
# iter <- 10 # number of iterations within EACH band, not overall
# nSite <- 10
# # centr <- FALSE; width <- 15
# centr <- TRUE; width <- 4
# absLat <- FALSE
# test <- bandit(dat, xy, width, iter, nSite, centr, absLat, output = 'full')

# TODO an option for equal-area latitudinal bands?

# cntr arg ignored if absLat = TRUE
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
