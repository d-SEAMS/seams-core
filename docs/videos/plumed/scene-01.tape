Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.32
# duration 16.10
# test docs_video_tapes
# test docs_engine_plumed
# test plumed_driver
# expect plumed hexagonal 128 six-rings 256 q6 0.625855
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 320ms
Type 'PLUMED_NUM_THREADS=1 python3 docs/howto/test_engines.py plumed . "$SEAMS_BUILD"'
Enter
Sleep 1s
Sleep 14780ms
