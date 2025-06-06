# ──────────────────────────────────────────────────────────────────────────────
# 1) builder stage: still do an editable install
# ──────────────────────────────────────────────────────────────────────────────
FROM mambaorg/micromamba:1.5.7 AS builder
ARG MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /opt/app

# Install OS build tools
USER root
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
        build-essential git && \
    rm -rf /var/lib/apt/lists/*
USER mambauser

# Create Conda environment
COPY environment.yml .
RUN micromamba create -y -n basebuddy -f environment.yml && \
    micromamba clean -afy

# Install pip-only dependency
RUN /opt/conda/envs/basebuddy/bin/pip install --no-cache-dir SigProfilerSimulator

# Copy source for editable install
# Explicitly copy data directory to ensure it's included
# COPY src/basebuddy/data /opt/app/src/basebuddy/data
COPY . .

# Ensure ownership so pip can write egg-info
USER root
RUN chown -R mambauser: /opt/app
USER mambauser

# Perform editable install (creates __editable__.basebuddy-*.pth)
RUN /opt/conda/envs/basebuddy/bin/pip install --no-cache-dir -e .

# ──────────────────────────────────────────────────────────────────────────────
# 2) runtime stage: copy both the Conda env and the source folder
# ──────────────────────────────────────────────────────────────────────────────
FROM debian:bookworm-slim
WORKDIR /work

# 1) Copy the Conda environment
COPY --from=builder /opt/conda /opt/conda

# 2) Copy the source tree (so the .pth reference still works)
COPY --from=builder /opt/app /opt/app

ENV PATH="/opt/conda/envs/basebuddy/bin:$PATH" \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

ENTRYPOINT ["basebuddy"]
CMD ["--help"]
