FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:ace9-main
#FROM python:3.9-slim-bullseye
RUN pip install latch==2.14.1
RUN mkdir /opt/latch

RUN apt-get update -y && \
    apt-get install -y curl unzip git

# Install Nextflow
RUN apt-get install -y default-jre-headless 
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/bin/ && \
    chmod 777 /usr/bin/nextflow 

# Install micromamba
ENV CONDA_DIR /opt/conda
ENV MAMBA_ROOT_PREFIX /opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
RUN apt-get update && apt-get install -y wget bzip2 \
    && wget -qO-  https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba \
    && touch /root/.bashrc \
    && ./bin/micromamba shell init -s bash -p /opt/conda  \
    && grep -v '[ -z "\$PS1" ] && return' /root/.bashrc  > /opt/conda/bashrc   # this line has been modified \
    && apt-get clean autoremove --yes \
    && rm -rf /var/lib/{apt,dpkg,cache,log}

SHELL ["bash", "-l" ,"-c"]

# Copy the original rnaseq-nf workflow 
COPY wf-nf /root/wf-nf

# Install conda dependencies
RUN micromamba create -f /root/wf-nf/conda.yml -y
ENV PATH=/opt/conda/envs/seqWell-nf/bin:$PATH

# Create workdir for results
RUN mkdir /root/work

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root