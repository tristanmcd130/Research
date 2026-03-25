FROM continuumio/miniconda3
WORKDIR /usr/local/app

RUN apt update && apt install -y graphviz
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority strict
RUN conda create -n sage sage python=3.11

SHELL ["conda", "run", "-n", "sage", "/bin/bash", "-c"]
RUN sage -pip install dot2tex

CMD ["conda", "activate", "sage"]