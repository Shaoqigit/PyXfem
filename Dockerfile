# use debian as base image
FROM debian:bullseye-slim

# install python3.11 and mini conda
RUN apt-get update && apt-get install -y python python3-pip

# install conda
RUN apt-get install -y wget
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh
RUN bash Miniconda3-py39_4.10.3-Linux-x86_64.sh -b
RUN rm Miniconda3-py39_4.10.3-Linux-x86_64.sh


# Set the environment variable
ENV PYTHONUNBUFFERED=1

# Create the working directory in the container
RUN mkdir /app

# Set the working directory
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app
# Copy the rest of your application's code
RUN ls -la
# install the dependencies with conda
RUN /root/miniconda3/bin/conda env create -f environment.yml

ENV PATH="/root/miniconda3/bin:${PATH}"

RUN conda init bash
# Use bash shell for RUN commands that need Conda environments
SHELL ["/bin/bash", "-c"]

# activate the conda environment
RUN echo "conda activate pyacoustix" >> ~/.bashrc

ENV BASH_ENV=~/.bashrc
# install other dependencies with pip
RUN /root/miniconda3/envs/pyacoustix/bin/pip install -r requirements.txt

# install package
RUN /root/miniconda3/envs/pyacoustix/bin/pip install .
