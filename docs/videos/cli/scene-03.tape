Output scene-03.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.40
# duration 13.27
# test docs_video_tapes
# expect cutoff 3.5
# expect nI 128
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 400ms
Type 'seams cn input/traj/ice1h.lammpstrj --cutoff 3.5 --type 1'
Enter
Sleep 1s
Sleep 11870ms
