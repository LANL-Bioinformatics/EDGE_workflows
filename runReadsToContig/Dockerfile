# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0 AS build

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda init bash \
    && . ~/.bashrc \
    && conda create --name readsToContig \
    && conda activate readsToContig 

RUN conda install -n readsToContig -c conda-forge r-base
RUN conda install -n readsToContig -c conda-forge python=3.11
RUN conda install -n readsToContig -c bioconda perl-json
RUN conda install -n readsToContig -c bioconda samclip=0.4.0
RUN conda install -n readsToContig -c bioconda bwa
RUN conda install -n readsToContig -c bioconda bowtie2
RUN conda install -n readsToContig -c bioconda samtools=1.6
RUN conda install -c conda-forge conda-pack


ADD bin/bam_to_fastq.pl /opt/conda/envs/readsToContig/bin
ADD bin/contig_stats.pl /opt/conda/envs/readsToContig/bin
ADD bin/ContigCoverageFold_plots_from_samPileup.pl /opt/conda/envs/readsToContig/bin
ADD bin/fastq_utility.pm /opt/conda/envs/readsToContig/bin
ADD bin/runReadsToContig.pl /opt/conda/envs/readsToContig/bin
ADD bin/tab2Json_for_dataTable.pl /opt/conda/envs/readsToContig/bin

RUN conda-pack -n readsToContig -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

RUN /venv/bin/conda-unpack

FROM debian:buster AS runtime
#RUN apk add --no-cache bash

COPY --from=build /venv /venv
ENV PERL5LIB=/venv/lib/perl5/core_perl
ENV PATH=/venv/bin:$PATH
    
SHELL ["/bin/bash", "-c"]
CMD /bin/bash