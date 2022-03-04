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

RUN mkdir -p data A_preprocessing B_annotation_GNPS C_annotation_SIRIUS

# Add scripts for running the container

ADD --chown=1000:1000 https://raw.githubusercontent.com/Coayala/MetaboTandem/main/scripts/A_preprocessing/1_data_import.ipynb /home/jovyan/work/A_preprocessing/1_data_import.ipynb
ADD --chown=1000:1000 https://raw.githubusercontent.com/Coayala/MetaboTandem/main/scripts/A_preprocessing/2_data_preprocessing.ipynb /home/jovyan/work/A_preprocessing/2_data_preprocessing.ipynb 
ADD --chown=1000:1000 https://raw.githubusercontent.com/Coayala/MetaboTandem/main/scripts/A_preprocessing/modified_writeMgfData.R /home/jovyan/work/A_preprocessing/modified_writeMgfData.R

ADD --chown=1000:1000 https://raw.githubusercontent.com/Coayala/MetaboTandem/main/scripts/B_annotation_GNPS/0_fixing_GNPS_output.ipynb /home/jovyan/work/B_annotation_GNPS/0_fixing_GNPS_output.ipynb

ADD --chown=1000:1000 https://raw.githubusercontent.com/Coayala/MetaboTandem/main/scripts/C_annotation_SIRIUS/1_ms2_into_groups.ipynb /home/jovyan/work/C_annotation_SIRIUS/1_ms2_into_groups.ipynb
ADD --chown=1000:1000 https://raw.githubusercontent.com/Coayala/MetaboTandem/main/scripts/C_annotation_SIRIUS/2_running_SIRIUS.ipynb /home/jovyan/work/C_annotation_SIRIUS/2_running_SIRIUS.ipynb
ADD --chown=1000:1000 https://raw.githubusercontent.com/Coayala/MetaboTandem/main/scripts/C_annotation_SIRIUS/3_fix_sirius_output.ipynb /home/jovyan/work/C_annotation_SIRIUS/3_fix_sirius_output.ipynb
ADD --chown=1000:1000 https://raw.githubusercontent.com/Coayala/MetaboTandem/main/scripts/C_annotation_SIRIUS/run_sirius_cli.sh /home/jovyan/work/C_annotation_SIRIUS/run_sirius_cli.sh

# Make executables scripts

RUN chmod +x /home/jovyan/work/C_annotation_SIRIUS/run_sirius_cli.sh
