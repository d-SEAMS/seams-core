Output scene-01.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 16.00
# test docs_video_tapes
# test lua_examples
# expect 2.14.0 cubic 8
# expect ffi ok
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/ice.lua .'
Enter
Sleep 1s
Type 'LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/luajit/?.so" luajit examples/lua/ice.lua .'
Enter
Sleep 1s
Type 'SEAMS_LIB=seams luajit tests/lua_ffi.lua "$SEAMS_BUILD/bindings/lua/dseams_ffi.lua"'
Enter
Sleep 13700ms
