# TODO - wrapping around antemeridian

# generate occurrences: 10 lat-long points in modern Australia
n <- 10
x <- seq(from = 140, to = 145, length.out = n)
y <- seq(from = -20, to = -25, length.out = n)
pts <- data.frame(x, y)

test_that('fails if points too sparse for any subsample', {
  expect_error(cookies(dat = pts, xy = 1:2, iter = 5, nSite = n+1, r = 200),
               'not enough close sites for any subsample')
  expect_error(clustr(dat = pts, xy = 1:2, iter = 5, nSite = n+1, distMax = 400),
               'insufficient points for any subsample')
  expect_error(bandit(dat = pts, xy = 1:2, bin = 30, nSite = n+1),
               'not enough close sites for any sample')
})

test_that('NA coordinates are removed/ignored', {
  pts[n, 1] <- NA # now there exist n-1 valid coordinates in data
  expect_error(cookies(dat = pts, xy = 1:2, iter = 5, nSite = n, r = 200),
               'not enough close sites for any subsample')
  expect_error(clustr(dat = pts, xy = 1:2, iter = 5, nSite = n, distMax = 400),
               'insufficient points for any subsample')
  expect_error(bandit(dat = pts, xy = 1:2, bin = 30, nSite = n),
               'not enough close sites for any sample')
})

test_that('multiple input classes ok', {
  expect_type(cookies(dat = pts, xy = 1:2,
                      iter = 5, nSite = n/2, r = 200),
              'list')
  expect_type(cookies(as.matrix(pts), xy = c('x', 'y'),
                      iter = 5, nSite = n/2, r = 200),
              'list')
  expect_type(clustr(dat = pts, xy = 1:2,
                     iter = 5, nSite = n/2, distMax = 400),
              'list')
  expect_type(clustr(as.matrix(pts), xy = c('x', 'y'),
                     iter = 5, nSite = n/2, distMax = 400),
              'list')
  expect_type(bandit(dat = pts, xy = 1:2,
                     iter = 5, nSite = n/2, bin = 10),
              'list')
  expect_type(bandit(as.matrix(pts), xy = c('x', 'y'),
                     iter = 5, nSite = n/2, bin = 10),
              'list')
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

test_that('projected coords ok', {
  expect_no_condition(cookies(dat = bivalves, xy = xyCell, crs = prj,
                              iter = 5, nSite = n, r = 1500))
  expect_no_condition(clustr( dat = bivalves, xy = xyCell, crs = prj,
                              iter = 5, nSite = n, distMax = 3000))
  expect_no_condition(bandit( dat = bivalves, xy = xyCell, crs = prj,
                              iter = 5, nSite = n, bin = 10, absLat = TRUE))
})

# TODO make bandit and clustr return error/warning!
test_that('coord projection flags pts outside lat-long bounds', {
  offPt <- xyFromCell(rPrj, 1)
  offPt[, 'y'] <- offPt[, 'y'] + 5000000 # shift north of north pole
  bivalves[1, xyCell] <- offPt # insert faulty coordinates into data
  expect_error(cookies(dat = bivalves, xy = xyCell, crs = prj,
                       iter = 5, nSite = n, r = 1500))
  # expect_warning()
  # clustr( dat = bivalves, xy = xyCell, crs = prj,
  #         iter = 5, nSite = n, distMax = 3000)
  #  bandit(dat = bivalves, xy = xyCell, crs = prj,
  #         iter = 10, nSite = 20, bin = 10)
})
