# mainly from https://pythonwise.blogspot.com/2015/04/docker-miniconda-perfect-match.html
FROM ubuntu:14.04
MAINTAINER MAGI team <magi_web@lbl.gov>
ENV DEBIAN_FRONTEND=noninteractive

# System packages 
# curl needed to install miniconda
# libXrender1 and libxext6 needed for rdkit
RUN apt-get update && apt-get install -y curl libXrender1 libxext6

# Install miniconda to /miniconda
ENV PATH=/miniconda/bin:${PATH}
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh && \
    bash Miniconda-latest-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda-latest-Linux-x86_64.sh
RUN conda update -y conda
COPY magi_env.yml /
RUN conda config --add channels conda-forge && \
    conda config --add channels rdkit && \
    conda env create -f magi_env.yml

# this is how you activate the environment
ENV PATH /miniconda/envs/magi/bin:$PATH

# Setup application
RUN mkdir magi
COPY workflow magi/workflow
COPY local_settings magi/local_settings
COPY tests magi/tests
COPY __init__.py magi/
COPY setup_docker.py magi/
RUN python /magi/setup_docker.py

ENTRYPOINT ["python", "/magi/workflow/magi_workflow.py"]

# commands:

# interactive shell
# docker run --entrypoint /bin/bash -it magi_test

# run a job
# docker run -v [LOCAL OUTPUTS DIR]:/magi/outputs -v [LOCAL INPUTS DIR]:/magi/inputs -t magi_test -f /magi/inputs/[FASTA INPUT] -c /magi/inputs/[COMPOUND INPUT]

# run a test
# docker run -v ~/Desktop/magi_outputs:/magi/outputs -v ~/Desktop/magi_inputs:/magi/inputs -t magi_test -f /magi/inputs/s_coelicolor_genes_fasta_smallset.faa -c magi/inputs/s_coelicolor_pactolus_data_smallset.csv