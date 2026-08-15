;; Full-API parity in Fennel: the same registered surface the Lua example
;; (example_lua/full_api/script.lua) drives, spelled idiomatically. Globals
;; come from vars.lua plus the trajectory path from config.yml.
(print "fennel_full: exercising the registered surface from Fennel")

(local scratch (.. (or (os.getenv "TMPDIR") "/tmp")
                   "/dseams_fennel_full_" (tostring (os.time)) "/"))

;; Readers -------------------------------------------------------------------
(local mw (readLammpsTrjO mwTrajectory targetFrame mwAtomType))
(assert (> mw.nop 0))
(local box (mw:box))
(print (string.format "read mw_nop=%d frame=%d box={%.2f, %.2f, %.2f}"
                      mw.nop mw.currentFrame (. box 1) (. box 2) (. box 3)))

;; Error propagation: a bad path raises a Lua error, caught by pcall ---------
(local (ok err) (pcall readXYZ "no/such/file.xyz"))
(assert (and (not ok) (not= err nil)))
(print "pcall_missing_file caught=true")

;; Neighbour lists -----------------------------------------------------------
(local nlMw (neighListO cutoffRadius mw mwAtomType))
(local nlMwIdx (neighbourListByIndex mw nlMw))
(assert (= (length nlMwIdx) mw.nop))
(print (string.format "neighbours rows=%d row1_len=%d"
                      (length nlMw) (length (. nlMwIdx 1))))

;; Rings: the incremental updater recomputes nothing on a repeat -------------
(local updater (RingUpdater.new maxDepth))
(local ringsFirst (updater:update nlMwIdx))
(local sourcesFirst (updater:lastRecomputedSources))
(local ringsSecond (updater:update nlMwIdx))
(assert (= (length ringsFirst) (length ringsSecond)))
(assert (= (updater:lastRecomputedSources) 0))
(print (string.format "rings n=%d first_sources=%d repeat_sources=%d"
                      (length ringsFirst) sourcesFirst
                      (updater:lastRecomputedSources)))

;; Bond-classification registry ----------------------------------------------
(local ruleNames (bondClassifierNames))
(var haveChill false)
(var haveChillPlus false)
(for [i 1 (length ruleNames)]
  (when (= (. ruleNames i) "CHILL") (set haveChill true))
  (when (= (. ruleNames i) "CHILL+") (set haveChillPlus true)))
(assert (and haveChill haveChillPlus))

(classifyBonds mw nlMw "CHILL")
(local chillTypes (getIceType mw nlMw scratch targetFrame false "chill.txt"))
(assert (= (length chillTypes) mw.nop))

(registerBondClassifier "fennel-staggered"
                        {:staggeredMax 1.0 :eclipsedMin 2.0 :eclipsedMax 3.0
                         :coordinationNumber 4})
(local mwCustom (readLammpsTrjO mwTrajectory targetFrame mwAtomType))
(classifyBonds mwCustom nlMw "fennel-staggered")
(print (string.format "bond_rules n=%d chill_by_name=true custom_registered=true"
                      (length ruleNames)))

;; Hydrogen bonds with explicit geometric thresholds -------------------------
;; Defaults are the water criterion (2.42 A, 30 deg); a 0.1 A cutoff sits
;; below any O-H separation and must produce an edge-free network
(local water (readLammpsTrjO trajectory targetFrame oxygenAtomType))
(local nlWater (neighListO cutoffRadius water oxygenAtomType))
(local hbDefault (getHbondNetwork trajectory water nlWater targetFrame 1))
(local hbNone (getHbondNetwork trajectory water nlWater targetFrame 1 0.1 30.0))
(var edges 0)
(var noneEdges 0)
(for [i 1 (length hbDefault)]
  (set edges (+ edges (- (length (. hbDefault i)) 1)))
  (set noneEdges (+ noneEdges (- (length (. hbNone i)) 1))))
(assert (and (> edges 0) (= noneEdges 0)))
(print (string.format "hbonds default_edges=%d tight_cutoff_edges=%d"
                      edges noneEdges))

;; The .con reader, when the readcon-core backend is compiled in -------------
(if (not= readCon nil)
    (do
      (local con (readCon conFile 2))
      (assert (= con.nop 4))
      (print (string.format "readcon nop=%d frame=%d" con.nop con.currentFrame)))
    (print "readcon backend absent; readCon not registered"))

;; Clean up the scratch output -----------------------------------------------
(os.remove (.. scratch "bop/chill.txt"))
(os.remove (.. scratch "bop"))
(os.remove scratch)

(print "fennel_full: done")
