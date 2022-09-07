context("correction")

test_that("left skewed", {
  N <- 1000
  
  set.seed(1)
  pvals <- runif(N, 0, 1)
  pvals <- sort(pvals)
  
  ind <- seq(1, 0.75 * N)
  pvals[ind] <- pvals[ind] / seq(5, 1, length = length(ind))
 
  pvals_c <- qq_correct(pvals)
  
  expect_true(qq_inflation(pvals) > 1)
  expect_equal(qq_inflation(pvals_c), 1, tolerance = .01)  
})
