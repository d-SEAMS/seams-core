Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 17.33
# test docs_video_tapes
# test docs_engine_gromacs
# test gromacs_mdrun
# expect gromacs hexagonal 128 six-rings 256 q6 0.625855
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'python3 docs/howto/test_engines.py gromacs . "$SEAMS_BUILD"'
Enter
Sleep 1s
Sleep 16030ms
