FROM jupyter/datascience-notebook:latest


# Install zip

USER root

RUN sudo apt update
RUN sudo apt install -y zip
RUN sudo apt-get install -y default-jre
RUN sudo apt install -y libnetcdf-*
RUN sudo apt-get install -y gcc-multilib g++-multilib

# Install R libraries required

RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')"
RUN R -e "install.packages('ncdf4', repos = 'https://cloud.r-project.org')"
RUN R -e "BiocManager::install('MSnbase')"
RUN R -e "BiocManager::install('mzR')"
RUN R -e "BiocManager::install('xcms')"

# Install Python libraries required

RUN pip install natsort

# Download SIRIUS

RUN wget https://bio.informatik.uni-jena.de/repository/dist-release-local/de/unijena/bioinf/ms/sirius/4.6.1/sirius-4.6.1-linux64-headless.zip
RUN unzip sirius-4.6.1-linux64-headless.zip
RUN echo "export SIRIUS=/home/jovyan/sirius/bin/sirius" >> /home/jovyan/.bashrc

# Create working directory

USER jovyan

WORKDIR /home/jovyan/work/

# Add scripts for running SIRIUS

