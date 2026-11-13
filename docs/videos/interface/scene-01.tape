Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 24.42
# test docs_video_tapes
# test docs_paper_interface_full
# expect interface full waters 2048
# expect interface full published hexagonal 840 interfacial 8 water 1131 clathrate 16 interClathrate 53
# expect interface full ovito hexagonal 807 interfacial 43 water 1140 clathrate 17 interClathrate 41
# expect interface full phi-c 0
# still interface/interface-full.lammpstrj
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'python3 docs/tutorials/paper_science.py interface-full "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/interface"'
Enter
Sleep 1s
Sleep 23120ms
