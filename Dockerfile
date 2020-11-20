################## BASE IMAGE #####################
FROM continuumio/miniconda3:4.7.12

################## METADATA #######################

LABEL base_image="continuumio/miniconda3"
LABEL version="4.7.12"
LABEL software="bam2cram"
LABEL software.version="1.0"
LABEL about.summary="Container image containing all requirements for save sam2cram conversion"
LABEL about.home="http://github.com/adigenova/bam2cram"
LABEL about.documentation="http://github.com/adigenova/bam2cram/README.md"
LABEL about.license_file="http://github.com/adigenova/bam2cram/LICENSE.txt"
LABEL about.license="GNU-3.0"


################## MAINTAINER ######################
MAINTAINER **digenovaa** <**digenovaa@fellows.iarc.fr**>


################## INSTALLATION ######################
COPY environment.yml /
RUN apt-get update && apt-get install -y procps && apt-get clean -y
RUN conda env create -n bam2cram -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/bam2cram/bin:$PATH
RUN conda env export --name bam2cram > bam2cram-v1.0.yml
