FROM condaforge/mambaforge
# NOT mambaforge-pypy3 (pandas is incompatible with pypy)

RUN apt-get update && apt-get install -y wget
# Download the bdbbind dataset
# RUN wget --no-clobber http://www.pdbbind.org.cn/download/PDBbind_v2020_refined.tar.gz
## update to the new download URL (around 10 times faster) from PDBbind website.
RUN wget --no-clobber https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v2020_refined.tar.gz
RUN tar -xvf PDBbind_v2020_refined.tar.gz

RUN mamba install -c conda-forge pandas

# The Docker image is generated on the fly, so we cannot pass any local files.
RUN wget --no-clobber -L https://raw.githubusercontent.com/PolusAI/mm-workflows/main/examples/scripts/generate_pdbbind_complex.py
RUN wget --no-clobber -L https://raw.githubusercontent.com/PolusAI/mm-workflows/main/examples/scripts/Dockerfile_pdbbind_refined
