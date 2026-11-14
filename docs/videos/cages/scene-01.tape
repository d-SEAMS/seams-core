Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.35
# duration 15.30
# test docs_video_tapes
# test docs_paper_cages_full
# expect cages full sI 16 48
# expect cages full sII 64 32
# expect cages full sH 24 16 8
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 350ms
Type 'python3 docs/tutorials/paper_science.py cages-full "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/cages"'
Enter
Sleep 1s
Sleep 13950ms
