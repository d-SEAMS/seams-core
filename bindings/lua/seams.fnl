;; Copyright (c) d-SEAMS developers
;; SPDX-License-Identifier: MIT
;;
;; Idiomatic Fennel over require("dseams").
;;
;; Errors are raised with error, the same convention as the Lua module.
;; A cloud is a table. Call (:close cloud) or let __gc free the frame.
;; Coordinates are a flat x, y, z sequence. The box is nine LAMMPS
;; numbers. Options are a table of keywords.

(local dseams (require :dseams))

(fn read [path ?opts]
  "Read one frame. ?opts.type is the species filter. A missing type keeps every atom when :all is true."
  (dseams.read path (or ?opts {})))

(fn system [xyz box ?types ?ids]
  "Build a frame from a flat coordinate table and a nine-number cell."
  (dseams.system xyz box ?types ?ids))

(fn chill+ [cloud ?opts]
  "CHILL+ labels. ?opts.rule is :published or :ovito. ?opts.graph is :cutoff, :mutual, or :seeded."
  (let [opts (or ?opts {})]
    (dseams.chill_plus cloud {:cutoff (or opts.cutoff 3.5)
                              :type (or opts.type 0)
                              :graph (or opts.graph :cutoff)
                              :rule (or opts.rule :published)})))

(fn chill [cloud ?opts]
  "Original CHILL labels. Options match chill+."
  (let [opts (or ?opts {})]
    (dseams.chill cloud {:cutoff (or opts.cutoff 3.5)
                         :type (or opts.type 0)
                         :graph (or opts.graph :cutoff)
                         :rule (or opts.rule :published)})))

(fn rings [cloud ?opts]
  "Primitive rings. ?opts.max-depth defaults to 6."
  (let [opts (or ?opts {})]
    (dseams.rings cloud {:cutoff (or opts.cutoff 3.5)
                         :type (or opts.type 0)
                         :graph (or opts.graph :cutoff)
                         :max_depth (or opts.max-depth opts.max_depth 6)})))

(fn neighbors [cloud ?opts]
  "Cutoff or mutual neighbour rows."
  (dseams.neighbors cloud (or ?opts {})))

(fn q4q6 [cloud ?opts]
  "Local q4 and q6 on the cutoff shell."
  (dseams.q4q6 cloud (or ?opts {})))

(fn cages [cloud ?opts]
  "Hexagonal and double-diamond membership. The default graph is seeded."
  (dseams.cages cloud (or ?opts {})))

{:read read
 :system system
 :chill+ chill+
 :chill chill
 :rings rings
 :neighbors neighbors
 :q4q6 q4q6
 :cages cages
 :version dseams.version
 :close (fn [cloud] (cloud:close))
 :core dseams.core
 :PUBLISHED dseams.PUBLISHED
 :OVITO dseams.OVITO}
