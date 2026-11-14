Output scene-02.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 12.00
# test docs_video_tapes
# expect nop 128 rings 256
# expect 6 256
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type '"$SEAMS" rings input/traj/ice1h.lammpstrj --type 1 --cutoff 3.5 --max-ring 6'
Enter
Sleep 1s
Sleep 10700ms
