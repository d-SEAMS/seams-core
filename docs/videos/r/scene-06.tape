Output scene-06.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.45
# duration 16.65
# test docs_video_tapes
# test docs_engine_openmm
# expect r example ok
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 450ms
Type 'Rscript docs/examples/minimal_r.R input/traj/genice_sI.lammpstrj'
Enter
Sleep 1s
Sleep 15200ms
