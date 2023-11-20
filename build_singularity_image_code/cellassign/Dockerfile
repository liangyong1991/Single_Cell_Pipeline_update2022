FROM condaforge/mambaforge:latest

COPY cellassign.yaml /opt/cellassign.yaml
COPY tensorflow-1.13.1-cp37-cp37m-manylinux1_x86_64.whl /opt/tensorflow-1.13.1-cp37-cp37m-manylinux1_x86_64.whl

RUN mamba env create -f /opt/cellassign.yaml --prefix /opt/conda/envs/cellassign
RUN /opt/conda/envs/cellassign/bin/python3.7 -m pip install /opt/tensorflow-1.13.1-cp37-cp37m-manylinux1_x86_64.whl

ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python
ENV PATH="/opt/conda/envs/cellassign/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/conda/envs/cellassign/lib:${LD_LIBRARY_PATH}"

RUN export PATH=/opt/conda/envs/cellassign/bin:${PATH}

RUN apt install tzdata
