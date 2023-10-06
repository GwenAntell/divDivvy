test_that('values return within basic tolerance limits', {
  smry <- sdSumry(bivalves, c('paleolng','paleolat'), 'genus', quotaQ = 0.8, quotaN = 100)
  expect_s3_class(smry, 'data.frame')
  expect_gt(smry[, 'SQSdiv'], 200)
  expect_gt(smry[, 'CRdiv'],  300)
  expect_gt(smry[, 'SQSupr95'], smry[, 'SQSlow95'])
  expect_gt(smry[, 'CRupr95'],  smry[, 'CRlow95'])
  expect_equal(smry[, 'nOcc'], 5096)
  expect_equal(smry[, 'meanPairDist'], 7277.259, tolerance = 0.01)
})

# project occurrences and retrieve cell centroids in new coordinate system
rWorld <- terra::rast()
prj <- 'EPSG:8857'
rPrj <- terra::project(rWorld, prj, res = 200000)
xyCartes <- c('paleolng','paleolat')
xyCell   <- c('centroidX','centroidY')
llOccs <- terra::vect(bivalves, geom = xyCartes, crs = 'epsg:4326')
prjOccs <- terra::project(llOccs, prj)
cellIds <- terra::cells(rPrj, prjOccs)[,'cell']
bivalves[, xyCell] <- terra::xyFromCell(rPrj, cellIds)

# insert faulty coordinates into data
offPt <- terra::xyFromCell(rPrj, 1)
offPt[, 'y'] <- offPt[, 'y'] + 5000000 # shift north of north pole
bivalves[1, xyCell] <- offPt

test_that('coord projection flags pts outside lat-long bounds', {
  expect_warning(expect_warning(expect_warning(
  expect_warning(expect_warning(expect_warning(
  expect_warning(sdSumry(bivalves, xyCell, 'genus') )
           #      'bounding box has potentially an invalid value range for longlat data')
    )))
  )))
})

