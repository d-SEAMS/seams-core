Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 14.00
# test docs_video_tapes
# expect nop 250 rings 511
# expect 3 37
# expect 4 273
# expect 5 34
# expect 6 167
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type '"$SEAMS" rings input/traj/exampleTraj.lammpstrj --type 2 --cutoff 3.5 --max-ring 6'
Enter
Sleep 1s
Sleep 12700ms
