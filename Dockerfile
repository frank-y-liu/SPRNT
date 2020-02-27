FROM ubuntu:bionic
MAINTAINER Frank Liu <https://github.com/frank-y-liu>
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -q \
    && apt-get install -qy build-essential wget curl libblas-dev git gfortran libgfortran-7-dev zlib1g zlib1g-dev bc \
    && rm -rf /var/lib/apt/lists/*

# build sprnt
RUN git clone https://github.com/frank-y-liu/SPRNT.git; \
    cd SPRNT; \
    ./configure; \
    make dep; \
    make; \
    make test; \
    make install; \    
    make clean
    
ENV PATH="/SPRNT/bin:${PATH}"

WORKDIR /data
VOLUME ["/data"]
