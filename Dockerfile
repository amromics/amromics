FROM mambaorg/micromamba:latest
#Metadata
LABEL \
	image.name="amromics" \
	image.version="0.1" \
	image.description="" \
	maintainer="amromics.org" 

# Set Environment Variables
#ENV PATH=$PATH:/opt/conda/bin

# Install amromics via micromamba
WORKDIR /tmp
RUN micromamba create -y -c conda-forge -c defaults --name amromics python=3.10 git && \
	eval "$(micromamba shell hook --shell bash)" && \
	micromamba activate amromics && \
	git clone --recursive https://github.com/amromics/amromics.git && \
	cd amromics && \
	micromamba install -y -c conda-forge -c bioconda -c anaconda -c etetoolkit -c rpetit3 -c defaults --file requirements.txt && \
	pip install . && \
	pip install panta && \
	micromamba list && \
	micromamba clean -afy

WORKDIR /tmp/amromics
RUN tar zxvf db.tar.gz

ENTRYPOINT ["/tmp/amromics/entrypoint.sh"]
#CMD ["amr-analysis.py","download_db"]
