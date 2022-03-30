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

RUN wget https://github.com/boecker-lab/sirius/releases/download/v4.9.12/sirius-4.9.12-linux64-headless.zip
RUN unzip sirius-4.9.12-linux64-headless.zip
RUN echo "export SIRIUS=/home/jovyan/sirius/bin/sirius" >> /home/jovyan/.bashrc

CMD ["start.sh"]

# Creating directory for scripts

USER jovyan

RUN git clone https://github.com/Coayala/MetaboTandem.git

# Create data and working directory

RUN mkdir /home/jovyan/data
WORKDIR /home/jovyan/work/

# Setting up container start-up

COPY start.sh /usr/local/bin
