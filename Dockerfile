FROM nfcore/base
MAINTAINER Andreas Wilm <wilma@gis.a-star.edu.sg>
LABEL authors="wilma@gis.a-star.edu.sg" \
    description="Docker image containing all requirements for the rpd-rnaseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/rpd-rnaseq-1.1dev/bin:$PATH
