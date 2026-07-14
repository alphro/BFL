test_that("package data and stan file exist", {
  expect_true(file.exists(system.file("extdata","phmrc_clean.csv", package="BFL")))
  expect_true(file.exists(system.file("stan","no_partial_labels.stan", package="BFL")))
})
