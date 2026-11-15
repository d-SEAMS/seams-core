# Copyright (c) d-SEAMS developers
library(dseams)

test_that("a missing path is an R error", {
  expect_error(read_frame("no-such-dump.lammpstrj"), regexp = ".")
})

test_that("an empty system stops with the handle text", {
  sys <- seams_system_new()
  expect_s3_class(sys, "dseams_system")
  expect_output(print(sys), "dseams system")
  expect_error(seams_chill(sys), regexp = ".")
})

test_that("summary methods print the frame", {
  root <- Sys.getenv("DSEAMS_ROOT", unset = "")
  tip <- file.path(root, "input", "traj", "exampleTraj.lammpstrj")
  skip_if(root == "" || !file.exists(tip), "DSEAMS_ROOT is not set")
  frame <- read_frame(tip, type = 2L)
  expect_output(summary(frame), "dseams frame")
  ice <- chill_plus(frame)
  expect_s3_class(ice, "dseams_ice")
  expect_output(print(ice), "CHILL\\+")
  expect_output(summary(ice), "water")
})

test_that("a ring census prints and summarises", {
  rings <- primitive_rings(read_frame(dseams_example("ice1h.lammpstrj"), type = -1L),
                           graph = "cutoff", max_depth = 6L)
  expect_s3_class(rings, "dseams_rings")
  expect_output(print(rings), "size")
  expect_output(summary(rings), "size")
})
