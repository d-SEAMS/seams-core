Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.35
# duration 16.00
# test docs_video_tapes
# expect nop 128 graph seeded signature 4:6,6:2 cages 64 atoms 128
# expect nop 128 graph seeded signature 6:7 cages 0 atoms 0
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 350ms
Type '"$SEAMS" cages input/traj/ice1h.lammpstrj --graph seeded --signature hc+ddc'
Enter
Sleep 1s
Sleep 14650ms
