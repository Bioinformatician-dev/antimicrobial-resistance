FROM python:3.12-slim

LABEL org.opencontainers.image.title="amr-consensus"
LABEL org.opencontainers.image.description="Multi-database consensus detection of antimicrobial resistance genes"

RUN apt-get update && \
    apt-get install -y --no-install-recommends ncbi-blast+ ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY pyproject.toml README.md ./
COPY amr_consensus/ ./amr_consensus/

RUN pip install --no-cache-dir .

WORKDIR /data
ENTRYPOINT ["amr-consensus"]
CMD ["--help"]
