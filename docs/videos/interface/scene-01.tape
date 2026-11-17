Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 12.33
# test docs_video_tapes
# test docs_paper_interface_full
# expect interface full input molecules 864 frames 1 file-frames 1 file input/traj/interface_f0.lammpstrj
# expect interface full paper molecules 864 frames 1
# expect interface full published hexagonal 352 interfacial 49 water 434 clathrate 3 interClathrate 25
# expect interface full stack HHHHHHHH....C.
# still interface/interface-full.lammpstrj
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'python3 docs/tutorials/paper_science.py interface-full "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/interface"'
Enter
Sleep 1s
Sleep 11030ms
