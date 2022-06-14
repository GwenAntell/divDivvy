# values of unknown env should be specified as NA, not e.g. 'unknown'
# dat can be df in which case xy and env should be col names
# or, omit dat arg and xy should be coords, env should be vector
# cutoff is inclusive, e.g. values >= cutoff will be assigned to class

#' @export
classRast <- function(grid, dat = NULL, xy, env, cutoff){
  if (! is.null(dat)){
    xy  <- dat[, xy]
    env <- dat[, env]
  }
  naEnv <- is.na(env)
  env <- env[!naEnv]
  xy  <-  xy[!naEnv,]

  # tally occurrences with each enviro within raster cells
  lvls <- unique(env)
  nLvl <- length(lvls)
  lyrL <- lapply(lvls, function(lvl){
    xyMat <- data.matrix( xy[env==lvl,] )
    terra::rasterize(xyMat, grid, fun = length)
  }  )
  lyrs <- terra::rast(lyrL)
  names(lyrs) <- lvls

  if (cutoff <= 0.5){ stop('cutoff must be greater than 0.5') }
  # set NA values to 0 so when summed across layers, result is not NA
  reclass <- matrix(c(NA, 0), ncol = 2)
  lyrsNoNA <- terra::classify(lyrs, reclass)
  prop <- lyrs / sum(lyrsNoNA) # convert to proportion
  cuts <- rbind(c(0, cutoff, 0), # binarise
                c(cutoff, 1, 1))
  # values = cutoff classified as above the cutoff (left-closed interval)
  mjtys <- terra::classify(prop, cuts,
                    include.lowest = TRUE,
                    right = FALSE)
  terra::values(grid) <- NA # ensure blank starting raster
  for (l in 1:nLvl){
    grid[ mjtys[[l]]==1 ] <- l
  }

  # change NA cells to indet. if they do have occurrence data
  # but no single environmental class exceeds cutoff
  anydat <- terra::rasterize( data.matrix(xy), grid)
  reclass <- cbind(1, nLvl+1)
  anydat <- terra::classify(anydat, reclass)
  grid <- terra::cover(grid, anydat)

  lvls <- c(lvls, 'indet.')
  attrs <- data.frame('num' = 1:nLvl)
  attrs$mainClass <- lvls
  terra::levels(grid) <- attrs
  return(grid)
}
