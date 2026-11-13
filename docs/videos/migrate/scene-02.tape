Output scene-02.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.35
# duration 20.30
# test docs_video_tapes
# test docs_migrate
# expect classes 89
# expect count 0
# expect largest 750
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 350ms
Type 'seams fingerprint input/traj/exampleTraj.lammpstrj --type 2'
Enter
Sleep 1s
Type 'seams pairs input/traj/exampleTraj.lammpstrj --site 1=cationHead,2=anion'
Enter
Sleep 1s
Type 'seams domains input/traj/exampleTraj.lammpstrj --site 1=cationHead,2=anion --subset polar'
Enter
Sleep 1s
Sleep 16950ms
