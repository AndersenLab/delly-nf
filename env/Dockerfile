FROM mambaorg/micromamba:1.5.0

RUN micromamba install -n base bcftools delly -c bioconda -y \
	&& micromamba clean -a -y

ARG MAMBA_DOCKERFILE_ACTIVATE=1
USER root

RUN apt-get --allow-releaseinfo-change update && \
    apt-get install -y procps && \
    apt-get clean && \
	rm -rf /var/lib/apt/lists/*


USER $MAMBA_USER
ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT source /usr/local/bin/_entrypoint.sh
