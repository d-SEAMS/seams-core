Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 13.17
# test docs_video_tapes
# test docs_migrate
# expect nop 250 hbonds 460
# expect cn 9.144
# expect cutoff 4.5
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'seams hbonds input/traj/exampleTraj.lammpstrj --type 2 --htype 1'
Enter
Sleep 1s
Type 'seams cn input/traj/exampleTraj.lammpstrj --pair-types 1,2 --cutoff 4.5'
Enter
Sleep 1s
Sleep 10870ms
