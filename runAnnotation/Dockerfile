# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:24.5.0-0 AS build

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

RUN conda init bash \
    && . ~/.bashrc \
    && conda create --name annotation \
    && conda activate annotation 

RUN conda install -n annotation -c bioconda perl-lwp-protocol-https
RUN conda install -n annotation -c conda-forge r-base
RUN conda install -n annotation -c bioconda prokka
RUN conda install -n annotation -c bioconda mummer
RUN conda install -n annotation -c bioconda blast=2.16
RUN conda install -n annotation -c bioconda perl-yaml

RUN conda install -c conda-forge conda-pack

#Custom implementation of RATT from EDGE
RUN wget https://ref-db.edgebioinformatics.org/EDGE/dev/edge_dev_thirdParty_softwares.tgz \
    && tar -xvzf edge_dev_thirdParty_softwares.tgz \
    && tar -xvzf thirdParty/RATT.tar.gz -C . \
    && mv RATT /opt/conda/envs/annotation/opt

#add scripts from this project to bin
ADD bin/* /opt/conda/envs/annotation/bin

#pack environment for runtime image
RUN conda-pack -n annotation -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar

RUN /venv/bin/conda-unpack

FROM debian:latest AS runtime

COPY --from=build /venv /venv

ENV PATH=/venv/bin:/venv/opt/RATT:$PATH

SHELL ["/bin/bash", "-c"]
CMD /bin/bash