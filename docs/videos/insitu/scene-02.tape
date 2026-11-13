Output scene-02.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.40
# duration 20.77
# test docs_video_tapes
# test docs_engine_openmm
# expect insitu hexagonal 128 six-rings 256
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 400ms
Type 'python3 docs/howto/test_engines.py insitu . "$SEAMS_BUILD"'
Enter
Sleep 1s
Sleep 19370ms
