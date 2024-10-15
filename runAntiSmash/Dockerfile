# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0 AS build

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda init bash \
    && . ~/.bashrc \
    && conda create --name sma \
    && conda activate sma 

# install dependencies for sma
RUN conda install -n sma -c bioconda antismash=7.1.0
RUN conda install -n sma -c bioconda bcbio-gff=0.7.0
RUN conda install -c conda-forge conda-pack

RUN conda-pack -n sma -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

RUN /venv/bin/conda-unpack

FROM debian:latest AS runtime

COPY --from=build /venv /venv
ENV PATH=/venv/bin:$PATH
    
SHELL ["/bin/bash", "-c"]
CMD /bin/bash