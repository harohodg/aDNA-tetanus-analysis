FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

# Make sure add-apt-repository is available
RUN apt update && \
	apt install -y software-properties-common

# Install GCC which makes it easier to install boost
RUN add-apt-repository ppa:ubuntu-toolchain-r/test && \
	apt update && \
	apt upgrade -y && \
	apt install --no-install-recommends -y gcc-10 \
	libgmp3-dev \
	git-all \
	python3 \
	python3-distro \
	wget

# install HTSlib prerequisites
RUN apt update && \
	apt install -y autoconf \
	cmake \
	automake \
	make \
	perl \
	zlib1g-dev \
	libbz2-dev \
	liblzma-dev \
	libcurl4-gnutls-dev \
	libssl-dev

# install HTSlib
RUN wget -O htslib-1.12.tar.bz2 https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 && \
	tar xvjf htslib-1.12.tar.bz2 && \
	cd htslib-1.12 && \
	autoreconf -i && \
	./configure && \
	make && \
	make install

# install Boost, then clean up intermediate files a little bit
RUN cd / && \
	wget -O boost_1_55_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz/download && \
	tar xzvf boost_1_55_0.tar.gz && \
	cd boost_1_55_0/ && \
	apt install -y --no-install-recommends build-essential \
	g++ \
	python-dev \
	autotools-dev \
	libicu-dev \
	libbz2-dev \
	libboost-all-dev && \
	apt autoremove --purge --yes && \
	apt clean && \
	rm -rf /vat/lib/apt/lists/*

# Install Octopus
RUN cd / && \
	git clone -b master https://github.com/luntergroup/octopus.git && \
	cd octopus && \
	./scripts/install.py

ENV PATH "$PATH:/octopus/bin/"
