;; Copyright (c) d-SEAMS developers
;; Eight-atom cubic cell through the Fennel module. Every atom is cubic.
(local seams (require :seams))
(local a (/ (* 2.75 4.0) (math.sqrt 3.0)))
(local fcc [[0 0 0] [0.5 0.5 0] [0.5 0 0.5] [0 0.5 0.5]])
(local basis [[0 0 0] [0.25 0.25 0.25]])
(local xyz [])
(each [_ site (ipairs fcc)]
  (each [_ shift (ipairs basis)]
    (table.insert xyz (* (+ (. site 1) (. shift 1)) a))
    (table.insert xyz (* (+ (. site 2) (. shift 2)) a))
    (table.insert xyz (* (+ (. site 3) (. shift 3)) a))))
(local cloud (seams.system xyz [0 0 0 a a a 0 0 0]))
(local labels ((. seams :chill+) cloud {:cutoff 3.5 :rule :published}))
(assert (= labels.counts.cubic 8) labels.counts.cubic)
(print (.. "fennel cubic " (tostring labels.counts.cubic)))
(seams.close cloud)
