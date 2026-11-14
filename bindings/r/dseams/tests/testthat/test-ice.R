ice1h <- dseams_example("ice1h.lammpstrj")
tip <- dseams_example("exampleTraj.lammpstrj")
nuc <- dseams_example("nucleation_f1.lammpstrj")

test_that("ice Ih is 128 hexagonal waters and 256 six-rings", {
  frame <- read_frame(ice1h, type = -1L)
  expect_s3_class(frame, "dseams_frame")
  expect_equal(nframes(ice1h), 1L)
  expect_equal(dim(positions(frame)), c(128L, 3L))
  ice <- chill_plus(frame, graph = "cutoff")
  expect_equal(sum(ice$name == "hexagonal"), 128L)
  expect_equal(attr(ice, "counts")[[2]], 128L)
  rings <- primitive_rings(frame, graph = "cutoff", max_depth = 6L)
  expect_equal(rings$n[rings$size == 6L], 256L)
  q <- q4_q6(frame)
  expect_equal(nrow(q), 128L)
  expect_true(all(is.finite(q$q6)))
  expect_length(symmetry(frame), 128L)
  expect_true(is.nan(jump_tau90(frame, frame, dt = 1)))
})

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
  # The default graph is the cutoff shell (seams_frame_steinhardt).
  # 0.389362 is the mutual four-nearest mean and must not be this value.
  got <- mean(q3$ql[is.finite(q3$ql)])
  expect_equal(got, 0.374897, tolerance = 1e-4)
  expect_false(abs(got - 0.389362) < 1e-4)
})

test_that("seeded cages under SEAMS_OFFLOAD match the host", {
  frame <- read_frame(nuc, type = 1L)
  saved_offload <- Sys.getenv("SEAMS_OFFLOAD", unset = NA)
  on.exit(if (is.na(saved_offload)) Sys.unsetenv("SEAMS_OFFLOAD") else
    Sys.setenv(SEAMS_OFFLOAD = saved_offload), add = TRUE)
  Sys.unsetenv("SEAMS_OFFLOAD")
  host <- cages(frame, graph = "seeded")
  Sys.setenv(SEAMS_OFFLOAD = "1")
  off <- cages(frame, graph = "seeded")
  Sys.unsetenv("SEAMS_OFFLOAD")
  expect_identical(off$cage, host$cage)
  for (key in c("n_ih", "n_ic", "n_both", "n_water")) {
    expect_equal(attr(off, key), attr(host, key))
  }
  n <- attr(host, "n_ih") + attr(host, "n_ic") + attr(host, "n_both") +
    attr(host, "n_water")
  expect_equal(n, nrow(host))
  expect_gt(attr(host, "n_both"), 0L)

  # SEAMS_BIN is a seams CLI from the same engine version.
  bin <- Sys.getenv("SEAMS_BIN", unset = "")
  skip_if(bin == "" || !file.exists(bin), "SEAMS_BIN is not set")
  line <- system2(bin, c("cages", nuc, "--graph", "seeded", "--cutoff", "3.5"),
                  stdout = TRUE)
  parts <- strsplit(line[[length(line)]], " +")[[1]]
  # nop N graph seeded hexagonal H cubic C both B water W
  expect_equal(as.integer(parts[[6]]), attr(host, "n_ih"))
  expect_equal(as.integer(parts[[8]]), attr(host, "n_ic"))
  expect_equal(as.integer(parts[[10]]), attr(host, "n_both"))
  expect_equal(as.integer(parts[[12]]), attr(host, "n_water"))
})

test_that("cubic ice is 4096 cubic and 4096 double-diamond", {
  # mW_cubic.lammpstrj (1.4 MB) is not shipped; DSEAMS_ROOT is a
  # seams-core checkout.
  root <- Sys.getenv("DSEAMS_ROOT", unset = "")
  cubic <- file.path(root, "input", "traj", "mW_cubic.lammpstrj")
  skip_if(root == "" || !file.exists(cubic), "DSEAMS_ROOT is not set")
  frame <- read_frame(cubic, type = 1L)
  ice <- chill_plus(frame)
  expect_equal(sum(ice$name == "cubic"), 4096L)
  found <- cages(frame)
  expect_equal(attr(found, "n_ic"), 4096L)
  expect_gte(nframes(cubic), 1L)
  expect_output(print(frame), "dseams frame")
  order <- q4_q6(frame)
  expect_equal(nrow(order), 4096L)
  tau <- jump_tau90(frame, frame, dt = 1)
  expect_true(is.na(tau) || is.finite(tau))
  expect_equal(nrow(steinhardt(frame, order_l = 6L, graph = "cutoff")), 4096L)
})
