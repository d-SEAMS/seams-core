Output scene-06.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.45
# duration 17.90
# test docs_video_tapes
# test docs_engine_openmm
# test docs_examples
# expect rust example ok
# expect q6 0.595
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 450ms
Type 'SEAMS_LIB_DIR="$SEAMS_BUILD_SRC" cargo run --offline --quiet --manifest-path bindings/rust/Cargo.toml --example minimal -- input/traj/genice_sI.lammpstrj'
Enter
Sleep 1s
Sleep 16450ms
