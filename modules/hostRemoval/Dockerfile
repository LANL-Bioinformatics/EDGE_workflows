# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0 AS build

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda init bash \
    && . ~/.bashrc \
    && conda create --name hostRemoval \
    && conda activate hostRemoval 

# install dependencies for hostRemoval
RUN conda install -n hostRemoval -c bioconda bwa=0.7.17
RUN conda install -n hostRemoval -c conda-forge r-base=4.4.0
RUN conda install -n hostRemoval -c bioconda seqkit
RUN conda install -c conda-forge conda-pack

ADD bin/host_reads_removal_by_mapping.pl /opt/conda/envs/hostRemoval/bin
ADD bin/removal_stats.pl /opt/conda/envs/hostRemoval/bin


RUN conda-pack -n hostRemoval -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

RUN /venv/bin/conda-unpack

FROM debian:latest AS runtime

COPY --from=build /venv /venv
RUN apt-get update && apt-get install procps -y && apt-get clean
ENV PERL5LIB=/venv/lib/perl5/core_perl
ENV PATH=/venv/bin:$PATH
    
SHELL ["/bin/bash", "-c"]
CMD /bin/bash