FROM continuumio/miniconda3:23.5.2-0 AS build

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda init bash \
    && . ~/.bashrc \
    && conda create --name assembly \
    && conda activate assembly 

# install dependencies
RUN conda install conda-libmamba-solver
RUN conda config --set solver libmamba
RUN conda install -n assembly -c conda-forge python=3.9
RUN conda install -n assembly -c bioconda samtools=1.17
RUN conda install -n assembly -c bioconda racon=1.5
RUN conda install -n assembly -c bioconda seqtk=1.3
RUN conda install -n assembly -c bioconda spades=3.15.5
RUN conda install -n assembly -c bioconda minimap2=2.26
RUN conda install -n assembly -c bioconda megahit=1.2.9
RUN conda install -n assembly -c bioconda idba=1.1.3
RUN conda install -n assembly -c bioconda unicycler=0.5.0
RUN wget https://github.com/ruanjue/wtdbg2/releases/download/v2.5/wtdbg-2.5_x64_linux.tgz \
    && tar -xvzf wtdbg-2.5_x64_linux.tgz \
    && cp wtdbg-2.5_x64_linux/* /opt/conda/envs/assembly/bin
RUN conda install -n assembly -c bioconda flye=2.9.2
RUN conda install -n assembly git
RUN conda install -c conda-forge conda-pack

ADD bin/extractLongReads.pl /opt/conda/envs/assembly/bin
ADD bin/renameFilterFasta.pl /opt/conda/envs/assembly/bin

RUN conda-pack -n assembly -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

RUN /venv/bin/conda-unpack

FROM debian:latest AS runtime

COPY --from=build /venv /venv

RUN /venv/bin/git clone --depth 1 https://gitlab.com/chienchi/long_read_assembly.git
RUN apt-get update && apt-get install procps -y && apt-get clean

ENV PATH="/venv/bin:/long_read_assembly:$PATH"
ENV PERL5LIB=/venv/lib/perl5/core_perl

SHELL ["/bin/bash", "-c"]
CMD /bin/bash