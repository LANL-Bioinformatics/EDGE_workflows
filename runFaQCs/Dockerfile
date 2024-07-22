# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0 AS build

ENV container=docker

# dependencies for faqcs
RUN apt-get update -y \
  && apt-get install -y zlib1g \
  && apt-get install -y jellyfish 

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda init bash \
  && . ~/.bashrc \
  && conda create --name runFaQCs \
  && conda activate runFaQCs 

# install faqcs
RUN conda install -n runFaQCs -c bioconda faqcs=2.10
#RUN conda install -n runFaQCs -c conda-forge jellyfish
#RUN conda install -n runFaQCs -c conda-forge zlib
RUN conda install -c conda-forge conda-pack

RUN conda-pack -n runFaQCs -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

RUN /venv/bin/conda-unpack

FROM debian:latest AS runtime

COPY --from=build /venv /venv
ENV PATH=/venv/bin:$PATH
    
SHELL ["/bin/bash", "-c"]
CMD /bin/bash

