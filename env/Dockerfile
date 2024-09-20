FROM mambaorg/micromamba:1.5.0

RUN micromamba install -n base bcftools delly -c bioconda -y \
	&& micromamba clean -a -y

ARG MAMBA_DOCKERFILE_ACTIVATE=1

USER $MAMBA_USER
ENV PATH="/opt/conda/bin:$PATH"
ENTRYPOINT source /usr/local/bin/_entrypoint.sh