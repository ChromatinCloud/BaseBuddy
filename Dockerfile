# Dockerfile

# --- Stage 1: Builder ---
# This stage uses micromamba to install all complex dependencies into a clean environment.
FROM mambaorg/micromamba:1.5.7 as builder

WORKDIR /opt/app

# Copy the project definition and environment files.
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml .
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml .
COPY --chown=$MAMBA_USER:$MAMBA_USER src ./src


# Create the Conda environment from the YAML file.
# This single command installs all conda and pip dependencies.
RUN micromamba create -y -n basebuddy -f environment.yml && \
    micromamba clean -afy

# --- Stage 2: Final Image ---
# This stage creates the final, smaller production image by copying
# only the necessary installed environment from the builder.
FROM debian:bookworm-slim

# Copy the entire Conda environment from the builder stage.
COPY --from=builder /opt/conda/envs/basebuddy /opt/conda/envs/basebuddy

# Configure the environment for the final image.
ENV PATH="/opt/conda/envs/basebuddy/bin:$PATH"
# Set the environment variable for BAMSurgeon to find the Picard JAR
# that was installed by the 'picard' conda package.
# Note: You may need to update the version number if the picard package changes.
ENV BAMSURGEON_PICARD_JAR="/opt/conda/envs/basebuddy/share/picard-3.1.1-0/picard.jar"

# Set a working directory for running the application.
WORKDIR /data

# Set the entrypoint, making the container act like the 'basebuddy' executable.
ENTRYPOINT ["basebuddy"]
# Provide a default command to run if no other is specified.
CMD ["--help"]
