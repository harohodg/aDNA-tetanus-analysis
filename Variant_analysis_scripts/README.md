# Discovery and analysis of TeNT variants

This subdirectory contains all the scripts needed to reproduce the variant analysis in Hodgkins et al.

## Prerequisites
See `INSTALL.md`. In short, you need a Linux-based operating system with `python` and `singularity` installed. `python` is needed for some processing scripts, while `singularity` is used to virtualize other software in reproducible containers.

N.B. Some other programs required, but virtually all Linux environments will have them; if your environment does not, you probably know what you're doing, so you're on your own ;)...

The `setup.sh` and `variants.sh` should work on 64-bit x86 Linux systems with the required software. I (@mjmansfi) have not tested any other software/hardware (e.g. WSL2, ARM) configurations.

## Input data
You will need to download the `read_alignments.tar.gz` archive, which is available at the following link:

[https://doi.org/10.6084/m9.figshare.23925768](https://doi.org/10.6084/m9.figshare.23925768)

The full variant analysis is wrapped in a single script, `variants.sh`. See `./scripts/variants.sh -h` for options.

## Expert options
The Singularity syntax in `variants.sh` is similar enough to Docker that it could be ported over. Older versions of these scripts supported both Docker and Singularity, but issues with mounting work directories, file paths, and permissions became too much of a hassle with Docker. It would not be difficult to copy and paste the same commands, replacing `singularity exec [IMAGE]` with Docker syntax like so:

`docker run -it [IMAGE] -v "${PWD}":/work -w /work [command]`

Note that whatever files the command requires need to be mounted in the working directory of the Docker container, so for some commands more than one mount flag (`-v [host dir]:[working dir in container]`, e.g., `-v "${PWD}"/../misc:/work/../misc`) is required.
