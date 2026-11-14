test_that("examples are listed and resolved", {
  files <- dseams_example()
  expect_true(all(c("ice1h.lammpstrj", "exampleTraj.lammpstrj") %in% files))
  expect_true(file.exists(dseams_example("ice1h.lammpstrj")))
  expect_error(dseams_example("missing.lammpstrj"), "no example")
})

test_that("bad input is an R error, not a crash", {
  expect_error(read_frame(tempfile(fileext = ".lammpstrj")))
  bad <- tempfile(fileext = ".lammpstrj")
  writeLines("not a dump", bad)
  on.exit(unlink(bad), add = TRUE)
  expect_error(nframes(bad))
  expect_equal(nrow(positions(read_frame(bad))), 0L)
})

test_that("print and positions agree on the atom count", {
  frame <- read_frame(dseams_example("ice1h.lammpstrj"), type = -1L)
  expect_output(print(frame), "nop 128")
  expect_equal(nrow(positions(frame)), 128L)
  expect_identical(colnames(positions(frame)), c("x", "y", "z"))
})

test_that("Steinhardt graphs and degrees return one row per atom", {
  frame <- read_frame(dseams_example("ice1h.lammpstrj"), type = -1L)
  for (l in c(3L, 4L, 6L)) {
    s <- steinhardt(frame, order_l = l)
    expect_s3_class(s, "data.frame")
    expect_equal(nrow(s), 128L)
    expect_identical(attr(s, "order_l"), l)
  }
  mutual <- steinhardt(frame, order_l = 6L, graph = "mutual")
  expect_equal(nrow(mutual), 128L)
})

test_that("the OVITO rule is selectable", {
  frame <- read_frame(dseams_example("exampleTraj.lammpstrj"), type = 2L)
  ovito <- chill_plus(frame, rule = "ovito-chill+")
  expect_equal(nrow(ovito), 250L)
  expect_error(chill_plus(frame, rule = "no-such-rule"))
})

test_that("threads can be set and read back", {
  old <- dseams_threads()
  on.exit(dseams_threads(old), add = TRUE)
  expect_gte(old, 1L)
  expect_equal(dseams_threads(1L), 1L)
  expect_error(dseams_threads(0L))
})

test_that("the engine version is the package version", {
  expect_identical(dseams_version(), as.character(utils::packageVersion("dseams")))
})
