# generate occurrence data
x  <- c(rep(1, 10), NA, NA)
y  <- c(rep(1, 5), NA, NA, 2:6)
sp <- c(rep(letters[1:3], 2),
        rep(letters[4:5], 3))
obs <- data.frame(x, y, sp)

# compare original and unique datasets:
# rows 4 and 5 removed as duplicates of 1 and 2, respectively
obs
uniqify(obs, taxVar = 3, xy = 1:2)

# using taxon identifications or other third variable is optional
uniqify(obs, xy = c('x', 'y'))

uniqify(obs, xy = c('x', 'y'), na.rm = FALSE)

# TODO test: what if na.rm = F but there are NA values?
# Works, but might be better to throw error - force user to deal with NAs
# Or, explain behavior and note contrast with subsampling functions
# (where any NA coordinates are removed/ignored)
