Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 11.08
# test docs_video_tapes
# test docs_paper_nucleation_full
# expect nucleation full input molecules 4096 frames 148 file-frames 148 file docs/.cache/figshare-11448702/nucleation.lammpstrj
# expect nucleation full paper molecules 4096 frames 148
# expect nucleation full frame 1 q6bar 0.26063
# expect nucleation full frame 148 q6bar 0.287256
# expect nucleation full lines 148
# still nucleation/nucleation-full.lammpstrj
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'python3 docs/tutorials/paper_science.py nucleation-full "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/nucleation"'
Enter
Sleep 1s
Sleep 9780ms
