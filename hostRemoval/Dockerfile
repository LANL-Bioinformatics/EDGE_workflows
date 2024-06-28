# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

# install dependencies for hostRemoval (bwa)
RUN conda install -c bioconda bwa=0.7.17

ADD host_reads_removal_by_mapping.pl /opt/conda/bin

CMD ["/bin/bash"]