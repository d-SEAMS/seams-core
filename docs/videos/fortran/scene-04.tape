Output scene-04.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.20
# duration 15.57
# test docs_video_tapes
# test docs_engine_openmm
# expect rings 256
# expect 6 256
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 200ms
Type 'seams rings input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1'
Enter
Sleep 1s
Sleep 14370ms
