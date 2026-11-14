Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 12.33
# test docs_video_tapes
# test docs_paper_brine_full
# expect brine full frame 0 hexagonal 93 water 152 front 18.50 stack HHHHH.......
# expect brine full frame 1 hexagonal 160 water 86 front 25.90 stack HHHHHHH.....
# expect brine full frame 2 hexagonal 218 water 38 front 37.00 stack HHHHHHHHHH..
# still brine/brine-full.lammpstrj
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'python3 docs/tutorials/paper_science.py brine-growth "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/brine"'
Enter
Sleep 1s
Sleep 11030ms
