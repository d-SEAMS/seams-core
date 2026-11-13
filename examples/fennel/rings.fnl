;; Copyright (c) d-SEAMS developers
;; Primitive rings on the oxygen network, the rings tutorial frame.
(local dseams (require :dseams))
(local root (or (. arg 1) "."))
(local cloud (dseams.read (.. root "/input/traj/exampleTraj.lammpstrj") {:type 2}))
(local rings (dseams.ring-network cloud 8))
(var line "tip rings")
(for [size 0 8]
  (set line (.. line " " (tostring (or (. rings.by_size size) 0)))))
(print line)
(cloud:close)
