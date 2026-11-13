# Read the GenIce sI cell, then the cutoff shell, q6, CHILL+, and rings.
path <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(nzchar(path))
library(dseams)
frame <- read_frame(path)
shell <- steinhardt(frame, order_l = 6L, cutoff = 3.5)
stopifnot(nrow(shell) == 46L)
stopifnot(mean(shell$ql) > 0.5)
ice <- chill_plus(frame)
stopifnot(sum(ice$name == "clathrate") == 46L)
rings <- primitive_rings(frame, graph = "cutoff")
six <- rings$n[rings$size == 6L]
five <- rings$n[rings$size == 5L]
stopifnot(five == 48L, six == 54L)
cat("r example ok\n")
