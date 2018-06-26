# use the ubuntu base image
FROM ubuntu:18.04

# maintainer
MAINTAINER Tobias Rausch rausch@embl.de

# install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    g++ \
    git \
    unzip \
    wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# install ATAC-seq pipeline using Bioconda
RUN cd /opt \
    && git clone https://github.com/tobiasrausch/ATACseq.git \
    && cd /opt/ATACseq/ \
    && make all

# add user
RUN groupadd -r -g 1000 ubuntu && useradd -r -g ubuntu -u 1000 -m ubuntu
USER ubuntu

# by default /bin/bash is executed
CMD ["/bin/bash"]
