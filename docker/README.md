# These are files for building a Docker image for MAGI. 

Basically, all you should need to do is pull master and then run `build.sh` from within this directory.

`build.sh` copies relevant files from the repo into this directory, then
builds the docker image, then deletes those files (to prevent this directory from
being unnecessarily large).

There are two linux-compiled BLAST binaries in here; do not delete these binaries!
They work with the linux distro that is built into the docker image.
