Output scene-02.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 16.00
# test docs_video_tapes
# expect nop 4096 graph seeded signature 4:6,6:2 cages 0 atoms 0
# expect nop 4096 graph seeded signature 6:7 cages 2048 atoms 4096
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type '"$SEAMS" cages input/traj/mW_cubic.lammpstrj --type 1 --signature hc+ddc'
Enter
Sleep 1s
Sleep 14700ms
