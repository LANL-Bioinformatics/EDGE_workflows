# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0 AS build

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda init bash \
  && . ~/.bashrc \
  && conda create --name runFaQCs \
  && conda activate runFaQCs 


# plotting and trimming software for long reads
RUN conda install -n runFaQCs -c bioconda nanoplot=1.43.0
RUN conda install -n runFaQCs -c bioconda porechop

#faqcs
RUN conda install -n runFaQCs -c bioconda faqcs=2.10

#dependencies for illumina_fastq_QC.pl
RUN conda install -n runFaQCs -c bioconda perl-parallel-forkmanager
RUN conda install -n runFaQCs -c bioconda perl-string-approx
RUN conda install -n runFaQCs -c conda-forge r-base

RUN conda install -c conda-forge conda-pack

ADD bin/* /opt/conda/envs/runFaQCs/bin

RUN conda-pack -n runFaQCs -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar


RUN /venv/bin/conda-unpack

FROM debian:latest AS runtime

COPY --from=build /venv /venv
RUN apt-get update && apt-get install procps -y && apt-get clean
ENV PATH=/venv/bin:$PATH
    
SHELL ["/bin/bash", "-c"]
CMD /bin/bash

