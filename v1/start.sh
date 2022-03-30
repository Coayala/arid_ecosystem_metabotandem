# Fixing permissions for /home/work/jovyan directory
# chown -R 1000:1000 /home/jovyan/work

# Copying scripts to /home/jovyan/work directory

cp -r /home/jovyan/MetaboTandem/scripts /home/jovyan/work/

# Starting Jupyter lab
jupyter lab
