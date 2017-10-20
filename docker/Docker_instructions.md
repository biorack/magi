# Using the MAGI Docker image
---
## Instructions:
1. Install [Docker](https://www.docker.com/)
1. Pull the image: `docker pull oerbilgin/magi`
1. Make at least one directory on your local machine to be mounted on the docker image, *e.g.* `~/Desktop/magi_inputs/` and `~/Desktop/magi_outputs/`
1. Move your input files into the directory you will use as inputs (*e.g.* `~/Desktop/magi_inputs/`)
1. Run the following docker command: `docker run -v [LOCAL OUTPUTS DIR]:/magi/outputs -v [LOCAL INPUTS DIR]:/magi/inputs -t magi_test -f /magi/inputs/[FASTA INPUT] -c /magi/inputs/[COMPOUND INPUT] [OTHER ARGS]`
    1. `[LOCAL OUTPUTS DIR]` and `[LOCAL INPUTS DIR]` are the local directories you made above
    1. `[FASTA INPUT]` and `[COMPOUND INPUT]` are the **file names** for your FASTA and compound input files, respectively
    1. On your local machine, `[FASTA INPUT]` and `[COMPOUND INPUT]` should be located inside `[LOCAL INPUTS DIR]`
    1. All MAGI files will be available on your local computer in `[LOCAL OUTPUTS DIR]`
    1. `[OTHER ARGS]` are all other arguments for MAGI, but please **do not use `-o` or `--outputs`**, since then you will never see your output files!
