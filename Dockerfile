FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for vcf2ethnicity pipeline"

COPY environment.yml /

RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/vcf2ethnicity-1.0/bin:/opt/fastngsadmix:$PATH

RUN apt-get -y update && apt-get -y install make wget g++ zlib1g zlib1g-dev  build-essential

COPY assets/fastngsadmix/source /opt/fastngsadmix

RUN cd /opt/fastngsadmix && make

