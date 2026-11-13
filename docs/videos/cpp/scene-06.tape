Output scene-06.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.45
# duration 21.23
# test docs_video_tapes
# test docs_engine_openmm
# expect nop 128 pairs 512 q6 0.625855 hexagonal 128 rings 256
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 450ms
Type 'g++ -Ibindings/cxx docs/videos/cpp/ice.cpp -L"$SEAMS_BUILD_SRC" -L"$SEAMS_LANGDIR" -lseams_lang -lseams -Wl,-rpath,"$SEAMS_BUILD_SRC" -Wl,-rpath,"$SEAMS_LANGDIR" -o "$VIDEO_OUT/ice-cpp" && "$VIDEO_OUT/ice-cpp" input/traj/ice1h.lammpstrj'
Enter
Sleep 1s
Sleep 19780ms
