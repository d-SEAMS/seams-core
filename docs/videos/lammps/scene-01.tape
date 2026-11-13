Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 19.83
# test docs_video_tapes
# test docs_tutorial_insitu
# expect lammps calls 5 stride 2 ice 320
# expect lammps hexagonal 0 six-rings 128 q6 0.628539
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'cd engines/lammps && "$SEAMS_BUILD/engines/run_hook" in.dseams'
Enter
Sleep 1s
Sleep 18530ms
