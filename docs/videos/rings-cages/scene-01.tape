Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 10.67
# test docs_video_tapes
# test docs_paper_rings_full
# expect rings full input molecules 4096 frames 1 file-frames 11 file input/traj/mW_cubic.lammpstrj
# expect rings full paper molecules 4096 frames 1 rings 8192
# expect rings full nop 4096 rings 8192
# expect rings full histogram 6 8192
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'python3 docs/tutorials/paper_science.py rings-full "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/rings"'
Enter
Sleep 1s
Sleep 9370ms
