# generate occurrence data
x  <- c(1, 1, 1, 3, NA)
y  <- c(2, 2, 2, 4, 5)
sp <- c('a', 'a', 'b', 'c', 'd')
obs <- data.frame(x, y, sp)

test_that('uniqify returns correct output', {
  u1 <- uniqify(obs, xy = 1:2, taxVar = 3)
  u2 <- uniqify(obs, xy = c('x', 'y'), taxVar = 'sp')
  u3 <- uniqify(obs, xy = 1:2)
  expect_identical(u1, obs[-c(2,5),] )
  expect_identical(u1, u2)
  expect_identical(u3, obs[-c(2,3,5),])
})

test_that('NAs retained if specified', {
  uNa <- uniqify(obs, xy = 1:2, taxVar = 3, na.rm = FALSE)
  expect_identical(uNa, obs[-2,])
})

test_that('taxa names as factors ok', {
  obs$sp <- as.factor(obs$sp)
  uFactr <- uniqify(obs, xy = 1:2, taxVar = 3)
  expect_identical(uFactr, obs[-c(2,5),] )
})
