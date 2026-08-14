;; Ring network on mW cubic ice, written in Fennel. Uses the same globals the
;; Lua examples see: trajectory, resCloud, cutoffRadius, oxygenAtomType, ...
(print "Fennel function evaluation environment")

(for [frame targetFrame finalFrame frameGap]
  (local cloud (readFrameOnlyOne trajectory frame resCloud oxygenAtomType
                                 false [0 0 0] [0 0 0]))
  (local nlist (neighborList cutoffRadius cloud oxygenAtomType))
  (local bonded (bondNetworkByIndex cloud nlist))
  (local rings (getPrimitiveRings bonded maxDepth))
  (print (string.format "fennel_ring_count %d" (length rings))))
