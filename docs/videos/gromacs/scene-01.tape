Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 22.00
# test docs_video_tapes
# test gromacs_mdrun
# expect gromacs frames 5 ice 320 atoms 320
# expect gromacs hexagonal 0 six-rings 128 q6 0.628519
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'engines/gromacs/run_mdrun.sh "$SEAMS_BUILD/engines/gmx-seams" "$SEAMS_BUILD_SRC"'
Enter
Sleep 1s
Sleep 20700ms
