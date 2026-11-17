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
Type 'gfortran -J "$VIDEO_OUT" -c bindings/fortran/src/seams_kinds.f90 -o "$VIDEO_OUT/seams_kinds.o" && gfortran -J "$VIDEO_OUT" -I "$VIDEO_OUT" -c bindings/fortran/src/seams_c.f90 -o "$VIDEO_OUT/seams_c.o" && gfortran -J "$VIDEO_OUT" -I "$VIDEO_OUT" -c bindings/fortran/src/seams_support.f90 -o "$VIDEO_OUT/seams_support.o" && gfortran -J "$VIDEO_OUT" -I "$VIDEO_OUT" -c bindings/fortran/src/seams.f90 -o "$VIDEO_OUT/seams.o" && gfortran -J "$VIDEO_OUT" -I "$VIDEO_OUT" -c bindings/fortran/src/seams_methods.f90 -o "$VIDEO_OUT/seams_methods.o" && gfortran -I "$VIDEO_OUT" docs/examples/minimal_fortran.f90 "$VIDEO_OUT/seams_kinds.o" "$VIDEO_OUT/seams_c.o" "$VIDEO_OUT/seams_support.o" "$VIDEO_OUT/seams.o" "$VIDEO_OUT/seams_methods.o" -o "$VIDEO_OUT/ice-f" -L"$SEAMS_BUILD_SRC" -lseams -Wl,-rpath,"$SEAMS_BUILD_SRC" && "$VIDEO_OUT/ice-f" input/traj/genice_sI.lammpstrj'
Enter
Sleep 1s
Sleep 19370ms
