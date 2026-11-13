Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 24.00
# test docs_video_tapes
# expect pixi run build
# expect meson setup build
# expect meson compile -C build
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'printf '"'"'%s\n'"'"' '"'"'pixi run build'"'"' '"'"'meson setup build'"'"' '"'"'meson compile -C build'"'"''
Enter
Sleep 1s
Sleep 22700ms
