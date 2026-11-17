Output scene-02.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.35
# duration 13.63
# test docs_video_tapes
# test docs_engine_lammps
# expect lammps hexagonal 128 six-rings 256 q6 0.625855
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 350ms
Type 'python3 docs/howto/test_engines.py lammps . "$SEAMS_BUILD"'
Enter
Sleep 1s
Sleep 12280ms
