Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 11.92
# test docs_video_tapes
# test docs_paper_nucleation_full
# expect nucleation full frame 0 q6bar 0.26828
# expect nucleation full frame 1 q6bar 0.287644
# expect nucleation full frame 2 q6bar 0.56248
# still nucleation/nucleation-full.lammpstrj
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'python3 docs/tutorials/paper_science.py nucleation-full "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/nucleation"'
Enter
Sleep 1s
Sleep 10620ms
