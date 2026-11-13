Output scene-02.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.25
# duration 10.20
# test docs_video_tapes
# test docs_engine_openmm
# expect nop 128
# expect frame 1
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 250ms
Type 'seams read input/traj/ice1h.lammpstrj'
Enter
Sleep 1s
Sleep 8950ms
