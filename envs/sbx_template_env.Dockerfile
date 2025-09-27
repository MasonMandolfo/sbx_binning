FROM condaforge/mambaforge:latest

# Setup
WORKDIR /home/sbx_binning_env

COPY envs/sbx_binning_env.yml ./

# Install environment
RUN conda env create --file sbx_binning_env.yml --name sbx_binning

ENV PATH="/opt/conda/envs/sbx_binning/bin/:${PATH}"

# "Activate" the environment
SHELL ["conda", "run", "-n", "sbx_binning", "/bin/bash", "-c"]

# Run
CMD "bash"