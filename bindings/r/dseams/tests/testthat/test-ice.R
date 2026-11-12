library(dseams)

root <- Sys.getenv("DSEAMS_ROOT", unset = "")
tip <- file.path(root, "input", "traj", "exampleTraj.lammpstrj")
cubic <- file.path(root, "input", "traj", "mW_cubic.lammpstrj")

skip_if(root == "" || !file.exists(tip), "DSEAMS_ROOT is not set")

test_that("CHILL+ on the mixed TIP4P frame matches the CLI", {
  frame <- read_frame(tip, type = 2L)
  expect_equal(nrow(positions(frame)), 250L)
  ice <- chill_plus(frame)
  expect_s3_class(ice, "data.frame")
  tab <- table(ice$name)
  expect_false("clathrate" %in% names(tab))
  expect_equal(as.integer(tab[["interClathrate"]]), 9L)
  expect_equal(as.integer(tab[["water"]]), 241L)
  expect_identical(ice$label, chill_plus(frame, graph = "cutoff")$label)
  mutual <- table(chill_plus(frame, graph = "mutual")$name)
  expect_equal(as.integer(mutual[["clathrate"]]), 1L)
  expect_equal(as.integer(mutual[["interClathrate"]]), 7L)
  expect_equal(as.integer(mutual[["water"]]), 242L)
  cages_df <- cages(frame)
  expect_equal(attr(cages_df, "n_water"), 250L)
  expect_equal(attr(cages_df, "n_ih"), 0L)
  rings <- primitive_rings(frame, max_depth = 8L)
  expect_equal(rings$n[rings$size == 6L], 224L)
  sig <- signature_cages(frame, "512")
  expect_true(all(sig >= 0L))
  expect_length(symmetry(frame), 250L)
  q3 <- steinhardt(frame, order_l = 3L)
  expect_equal(mean(q3$ql[is.finite(q3$ql)]), 0.389362, tolerance = 1e-4)
})

test_that("seeded cages under SEAMS_OFFLOAD match the host and the CLI", {
  nuc <- file.path(root, "input", "traj", "nucleation_f1.lammpstrj")
  skip_if(!file.exists(nuc))
  frame <- read_frame(nuc, type = 1L)
  Sys.unsetenv("SEAMS_OFFLOAD")
  host <- cages(frame, graph = "seeded")
  Sys.setenv(SEAMS_OFFLOAD = "1")
  off <- cages(frame, graph = "seeded")
  Sys.unsetenv("SEAMS_OFFLOAD")
  expect_identical(off$cage, host$cage)
  expect_equal(attr(off, "n_ih"), attr(host, "n_ih"))
  expect_equal(attr(off, "n_ic"), attr(host, "n_ic"))
  expect_equal(attr(off, "n_both"), attr(host, "n_both"))
  expect_equal(attr(off, "n_water"), attr(host, "n_water"))
  n <- attr(host, "n_ih") + attr(host, "n_ic") + attr(host, "n_both") +
    attr(host, "n_water")
  expect_equal(n, nrow(host))
  expect_gt(attr(host, "n_both"), 0L)
  bin <- Sys.getenv("SEAMS_BIN", unset = "")
  skip_if(bin == "" || !file.exists(bin))
  # R prepends its own lib directories to LD_LIBRARY_PATH. The engine is
  # linked against the toolchain libc, so those prefixes are removed for
  # this one call.
  saved_ld <- Sys.getenv("LD_LIBRARY_PATH", unset = "")
  parts <- strsplit(saved_ld, ":", fixed = TRUE)[[1]]
  parts <- parts[!grepl("^/usr/lib", parts)]
  Sys.setenv(LD_LIBRARY_PATH = paste(parts, collapse = ":"))
  line <- system2(bin, c("cages", nuc, "--graph", "seeded", "--cutoff", "3.5"),
                  stdout = TRUE)
  Sys.setenv(LD_LIBRARY_PATH = saved_ld)
  parts <- strsplit(line[[length(line)]], " +")[[1]]
  # nop N graph seeded hexagonal H cubic C both B water W
  expect_equal(as.integer(parts[[6]]), attr(host, "n_ih"))
  expect_equal(as.integer(parts[[8]]), attr(host, "n_ic"))
  expect_equal(as.integer(parts[[10]]), attr(host, "n_both"))
  expect_equal(as.integer(parts[[12]]), attr(host, "n_water"))
})

test_that("cubic ice is 4096 cubic and 4096 double-diamond", {
  skip_if(!file.exists(cubic))
  frame <- read_frame(cubic, type = 1L)
  ice <- chill_plus(frame)
  expect_equal(sum(ice$name == "cubic"), 4096L)
  found <- cages(frame)
  expect_equal(attr(found, "n_ic"), 4096L)
})
