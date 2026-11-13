Output scene-02.mp4
Set Shell bash
Set Width 1920
Set Height 1080
Set Framerate 30
Set Padding 0
# lead 0.22
# duration 11.00
# test docs_video_tapes
# test docs_paper_brine_full
# expect brine temperature 240
Hide
Type 'PS1='"'"'$ '"'"''
Enter
Show
Sleep 220ms
Type 'python3 -c "import pathlib; text=pathlib.Path('"'"'docs/tutorials/paper_science.py'"'"').read_text(); assert '"'"'240 * unit.kelvin'"'"' in text; print('"'"'brine temperature 240'"'"')"'
Enter
Sleep 1s
Sleep 9780ms
