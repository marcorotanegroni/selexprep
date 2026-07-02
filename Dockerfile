# syntax=docker/dockerfile:1
#
# selexprep container image.
#
# Build:   docker build -t selexprep:0.2.0 .
# Run:     docker run --rm selexprep:0.2.0 inspect SRR1234567
#          docker run --rm -v "$PWD":/data selexprep:0.2.0 run /data/accessions.tsv --outdir /data/out
#
# The image installs selexprep from the local source tree (this checkout), so it
# is buildable before the PyPI release. cutadapt arrives as a declared pip
# dependency and is invoked as a subprocess at runtime — no apt package needed
# (its manylinux wheel runs on the debian-slim base).

FROM python:3.12-slim

LABEL org.opencontainers.image.title="selexprep" \
      org.opencontainers.image.description="Accession-first preprocessing for public HT-SELEX with primer auto-inference" \
      org.opencontainers.image.source="https://github.com/marcorotanegroni/selexprep" \
      org.opencontainers.image.licenses="MIT" \
      org.opencontainers.image.version="0.2.0"

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /build

# Copy only the files the build backend (hatchling) needs to produce the wheel.
# pyproject references README.md (readme) and LICENSE (license-files); the
# package code and bundled catalog CSV live under src/.
COPY pyproject.toml README.md LICENSE ./
COPY src ./src

RUN python -m pip install --upgrade pip \
    && python -m pip install . \
    && rm -rf /build

# Run as an unprivileged user; default working dir is a writable home.
RUN useradd --create-home --uid 1000 selex
USER selex
WORKDIR /home/selex

ENTRYPOINT ["selexprep"]
CMD ["--help"]
