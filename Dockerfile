FROM ubuntu:16.04

# install dependencies from pip3

RUN apt update && \
    apt install -y python3 ncbi-blast+ && \
    apt install -y python-biopython \
                   python3-pip \
                   wget \
                   unzip

# Install dependencies from conda 
RUN cd /usr/local/ && \
    wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda && \
    rm Miniconda3-latest-Linux-x86_64.sh && \
    ln -s /usr/local/miniconda/bin/conda /usr/local/bin/ && \
    conda init bash && \
    /bin/bash -c "source /root/.bashrc" && \
    conda install -c bioconda seqtk bowtie2 krakenuniq kallisto gmap snap-aligner openssl=1.0 samtools=1.7 bedtools bwa mafft bcftools tabix && \
    conda clean -afy
# Install Picard 

# Install PrimerClip 
RUN wget -qO- https://get.haskellstack.org/ | sh
#RUN /usr/local/bin/stack build 
RUN wget https://github.com/swiftbiosciences/primerclip/archive/deltest.zip && unzip deltest.zip && cd /primerclip-deltest/ && /usr/local/bin/stack build&&  /usr/local/bin/stack install && cd .. 

###########
# ANNOVAR #
###########


# http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz




##########
# JAVA 8 #
##########

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN export JAVA_HOME



CMD ["/bin/bash"]