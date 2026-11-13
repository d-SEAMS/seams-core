;; Copyright (c) d-SEAMS developers
;; CHILL+ on the basal ice-water interface, the cutoff list, published rule.
(local dseams (require :dseams))
(local root (or (. arg 1) "."))
(local cloud (dseams.read (.. root "/input/traj/interface_f0.lammpstrj") {:type 1}))
(local labels (dseams.chill-plus cloud {:cutoff 3.5 :type 1 :graph "cutoff" :rule dseams.PUBLISHED}))
(local order [:cubic :hexagonal :water :interfacial :clathrate :interClathrate :unclassified :reCubic :reHex])
(var line "interface chill")
(each [_ name (ipairs order)]
  (set line (.. line " " (tostring (. labels.counts name)))))
(print line)
(cloud:close)
