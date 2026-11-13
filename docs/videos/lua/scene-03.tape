Output scene-03.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.35
# duration 18.00
# test docs_video_tapes
# test lua_examples
# expect hexagonal 128
# expect q6bar 0.562448843702
# expect cages 128 0 0 0
# expect knn 512
# expect sheet 256
# expect sig 2 40
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 350ms
Type 'LUA_PATH=bindings/lua/?.lua LUA_CPATH="$SEAMS_BUILD/bindings/lua/lua54/?.so" lua5.4 docs/examples/lua_analyses.lua'
Enter
Sleep 1s
Sleep 16650ms
