# Workflow tools and PETRA itself. Databases, including the NCBI taxonomy used
# by taxonkit, are not part of the image: `petra run --use-apptainer` mounts the
# database directory built by `petra build-db`.
FROM mambaorg/micromamba:1.5-jammy

RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.12 \
    pip \
    diamond \
    taxonkit \
    cd-hit \
    minimap2 \
    samtools \
    bedtools \
    blast \
    pyrodigal \
    infernal \
    trnascan-se \
    && micromamba clean --all --yes

# Whokaryote pins Python 3.8, below what PETRA needs, so it gets its own
# environment. The wrapper puts that environment first on PATH so Whokaryote
# finds its own Tiara and Prodigal.
RUN micromamba create -y -n whokaryote -c conda-forge -c bioconda whokaryote \
    && micromamba clean --all --yes
USER root
RUN printf '#!/bin/sh\nPATH="/opt/conda/envs/whokaryote/bin:$PATH" exec /opt/conda/envs/whokaryote/bin/whokaryote.py "$@"\n' \
    > /usr/local/bin/whokaryote.py \
    && chmod 755 /usr/local/bin/whokaryote.py
USER $MAMBA_USER
ARG MAMBA_DOCKERFILE_ACTIVATE=1
COPY --chown=$MAMBA_USER:$MAMBA_USER . /tmp/petra
RUN pip install --no-cache-dir /tmp/petra && rm -rf /tmp/petra

ENV PATH="/opt/conda/bin:$PATH"
