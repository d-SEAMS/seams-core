Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 15.25
# test docs_video_tapes
# expect pixi run build
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'printf '"'"'%s\n'"'"' '"'"'pixi run build'"'"' '"'"'meson setup build'"'"' '"'"'meson compile -C build'"'"''
Enter
Sleep 1s
Sleep 13950ms
