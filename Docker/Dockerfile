################################################
# Dockerfile to build BioToolBox in a container
# Based on Ubuntu
################################################


FROM ubuntu:18.04

# Update aptitude and install some packages
# also install some of the known standard Perl dependencies 
# from ubuntu repos to make things a little easier
RUN apt-get update && apt-get -y install \
    build-essential \
    perl \
    perl-base \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    openssl \
    curl \
    libcurl4 \
    libcurl4-openssl-dev \
    cpanminus \
    libmodule-build-perl \
    libparallel-forkmanager-perl \
    libset-intervaltree-perl \
    libstatistics-descriptive-perl \
    libstatistics-lite-perl \
    libarchive-zip-perl \
    libtest-most-perl \
    libio-string-perl \
    libdata-types-perl \
    libdata-swap-perl \
    libfile-which-perl \
    libconfig-simple-perl \
    libtemplate-tiny-perl


# Set up environment
ENV OPT_SRC /opt/src


# Download files
WORKDIR $OPT_SRC
RUN curl -o htslib-1.9.tar.bz2 \
         -L https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    curl -o libBigWig-0.4.4.tar.gz \
         -L https://github.com/dpryan79/libBigWig/archive/0.4.4.tar.gz && \
    tar -xf htslib-1.9.tar.bz2 && \
    tar -xf libBigWig-0.4.4.tar.gz

# Install HTS library
WORKDIR $OPT_SRC/htslib-1.9
RUN make && make install

# Install libBigWig
WORKDIR $OPT_SRC/libBigWig-0.4.4
RUN make && make install


# Install Perl modules
# We are explicitly avoiding the current BioPerl version as well as the Apt version.
# Also turning off testing to speed up and simplify.
WORKDIR $OPT_SRC
RUN cpanm --notest \
    https://cpan.metacpan.org/authors/id/C/CJ/CJFIELDS/BioPerl-1.007002.tar.gz \
    https://cpan.metacpan.org/authors/id/A/AV/AVULLO/Bio-DB-HTS-3.01.tar.gz \
    https://github.com/Ensembl/Bio-DB-Big/archive/master.zip \
    Set::IntSpan::Fast::XS \
    https://cpan.metacpan.org/authors/id/T/TJ/TJPARNELL/Bio-ToolBox-1.67.tar.gz


# Install UCSC utilities
RUN curl -O -L http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig && \
    curl -O -L http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToWig && \
    curl -O -L http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed && \
    chmod +x wigToBigWig bigWigToWig bedToBigBed && \
    mv wigToBigWig bigWigToWig bedToBigBed /usr/local/bin/


# cleanup
RUN rm -rf $OPT_SRC /root/.cpanm
WORKDIR /data

