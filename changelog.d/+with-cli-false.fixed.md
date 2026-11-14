`-Dwith_cli=false` configures again; `seams_lang_fixtures` referenced the CLI target unconditionally, so `pip install .` failed.
