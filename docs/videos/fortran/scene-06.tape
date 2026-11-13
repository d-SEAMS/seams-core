Output scene-06.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.45
# duration 20.82
# test docs_video_tapes
# test docs_engine_openmm
# expect fortran example ok
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 450ms
Type 'gfortran -J "$VIDEO_OUT" -c bindings/fortran/seams.f90 -o "$VIDEO_OUT/seams.o" && gfortran -I "$VIDEO_OUT" docs/examples/minimal_fortran.f90 "$VIDEO_OUT/seams.o" -o "$VIDEO_OUT/ice-f" -L"$SEAMS_BUILD_SRC" -lseams -Wl,-rpath,"$SEAMS_BUILD_SRC" && "$VIDEO_OUT/ice-f" input/traj/genice_sI.lammpstrj'
Enter
Sleep 1s
Sleep 19370ms
