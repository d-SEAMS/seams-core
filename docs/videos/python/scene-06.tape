Output scene-06.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.45
# duration 19.57
# test docs_video_tapes
# test docs_engine_openmm
# expect nop 128 nI 128 q6 0.625855 hexagonal 128 rings 256
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 450ms
Type 'python3 docs/videos/python/ice.py input/traj/ice1h.lammpstrj "$SEAMS_BUILD_SRC"'
Enter
Sleep 1s
Sleep 18120ms
