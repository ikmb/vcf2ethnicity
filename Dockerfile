FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for vcf2ethnicity pipeline"

ADD . .

RUN mkdir -p /opt/ && cp -R /assets/fastngsadmix /opt/

RUN apt-get -y update && apt-get -y install make wget g++ zlib1g zlib1g-dev  build-essential ruby-full ruby-dev

RUN cd /opt/fastngsadmix && make

RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/vcf2ethnicity-1.0/bin:/opt/fastngsadmix:$PATH

