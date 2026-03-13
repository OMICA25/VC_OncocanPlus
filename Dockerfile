FROM continuumio/miniconda3:latest

# Copy your environment file
COPY environment.yml /tmp/environment.yml

# Create environment defined in environment.yml
RUN conda env create -f /tmp/environment.yml

# Activate by default when container starts
SHELL ["/bin/bash", "-c"]
RUN echo "source activate ngs_env" >> ~/.bashrc

# Ensure tools from the environment are in PATH
ENV PATH /opt/conda/envs/ngs_env/bin:$PATH

# Clean conda caches to reduce image size
RUN conda clean -afy
