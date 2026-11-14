Output scene-02.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.22
# duration 8.92
# test docs_video_tapes
# test docs_paper_brine_full
# expect brine waters 864 ions 9
# expect brine full input molecules 873 frames 1 file-frames 1 file input/traj/nacl_front.lammpstrj
# expect brine full paper waters 864 ions 9 frames 1
# expect brine full waters 864 ions 9 hexagonal 352 water 434
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 220ms
Type 'python3 -c "import pathlib; text=pathlib.Path('"'"'docs/tutorials/paper_science.py'"'"').read_text(); assert '"'"'brine full paper waters 864 ions 9 frames 1'"'"' in text; print('"'"'brine waters 864 ions 9'"'"')"'
Enter
Type 'python3 docs/tutorials/paper_science.py brine-full "$SEAMS" "$SEAMS_LIBDIR" "$VIDEO_OUT/brine"'
Enter
Sleep 1s
Sleep 7700ms
