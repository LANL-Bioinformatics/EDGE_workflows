# syntax=docker/dockerfile:1
FROM continuumio/miniconda3:23.5.2-0

ENV container=docker

# add conda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda

# install dependencies
RUN conda install -c bioconda phame=1.0.2
RUN conda install -c conda-forge libgcc-ng=12
RUN conda install -c bioconda minimap2 --solver classic
RUN conda install -c bioconda samclip=0.4.0
RUN conda install -c conda-forge pysam
RUN conda install bcftools "conda-forge::blas=*=openblas"
RUN conda install -c conda-forge r-base=4.4.0

ADD runReadsToGenome.pl /opt/conda/bin
ADD fastq_utility.pm /opt/conda/bin
ADD vcfutils.pl /opt/conda/bin
ADD align_trim.py /opt/conda/bin

CMD ["/bin/bash"]