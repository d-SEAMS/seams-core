Output scene-04.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.40
# duration 14.00
# test docs_video_tapes
# test lua_examples
# expect interface chill 1 352 434 49 3 25 0 0 0
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 400ms
Type 'lua5.4 tests/lua_fennel_compile.lua "$FENNEL_LUA" examples/fennel/interface.fnl > "$VIDEO_OUT/interface.lua" && LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 "$VIDEO_OUT/interface.lua" .'
Enter
Sleep 1s
Sleep 12600ms
