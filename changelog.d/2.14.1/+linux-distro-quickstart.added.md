The release upload now checks the linux-x86_64 tarball in fresh Ubuntu
22.04, Rocky Linux 9, and Debian 12 containers. Each check compiles
`docs/examples/quickstart.c` with the `pkg-config` flags for `seams` and
runs it on an eight-molecule ice frame.
