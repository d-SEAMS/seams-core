Output scene-02.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.30
# duration 18.00
# test docs_video_tapes
# test lua_examples
# expect dseams.read nop=250
# expect clathrate 46
# expect neighbors rows=46 nnz=184
# expect q6 0.595133082718
# expect rings 3 37 4 273 5 34 6 167 7 0 8 0 fraction 1.000000
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 300ms
Type 'LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/read.lua .'
Enter
Sleep 1s
Type 'LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/chill_plus.lua .'
Enter
Sleep 1s
Type 'LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/neighbors.lua .'
Enter
Sleep 1s
Type 'LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/steinhardt.lua .'
Enter
Sleep 1s
Type 'LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 examples/lua/rings.lua .'
Enter
Sleep 12700ms
