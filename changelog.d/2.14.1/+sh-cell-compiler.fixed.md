`docs_sh_cell` compiles its header snippet with the C compiler the build was configured with, and without the suite's library path, so the loader opens the libc that compiler linked.
