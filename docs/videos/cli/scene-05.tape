Output scene-05.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.35
# duration 11.13
# test docs_video_tapes
# test docs_engine_openmm
# expect q6 0.625855
# expect hexagonal 128
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 350ms
Type 'seams steinhardt input/traj/ice1h.lammpstrj --cutoff 3.5 --l 6 --type 1'
Enter
Sleep 1s
Type 'seams chill-plus input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1'
Enter
Sleep 1s
Sleep 8780ms
