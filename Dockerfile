# Essential
FROM ubuntu:20.04
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -yq apt-utils dialog
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get install -yq build-essential software-properties-common
RUN apt-get update

# acados dependencies
RUN apt-get install -yq git cmake make
RUN apt-get install -yq bc libblas-dev liblapack-dev

# Deep bug somwhere when installing matplotlib on ubuntu20
# https://stackoverflow.com/questions/25674612/ubuntu-14-04-pip-cannot-upgrade-matplotllib
RUN apt-get install -yq libfreetype6-dev libxft-dev

# python3.8
RUN add-apt-repository -y ppa:deadsnakes/ppa
RUN apt-get -y update
RUN apt-get install -y python3.8 python3.8-tk python3.8-dev
RUN apt-get -y install python3-pip
RUN /usr/bin/python3.8 -m pip install --upgrade pip
RUN pip3 install --upgrade pip


# Install python dependencies
RUN mkdir -p /root/src/
WORKDIR /root/src
COPY . .

#ACADOS
ENV LD_LIBRARY_PATH=/root/src/acados/lib
ENV ACADOS_SOURCE_DIR=/root/src/acados
SHELL ["/bin/bash", "-c"]

RUN pip3 install virtualenv
RUN virtualenv acadosenv --python=python3.8
RUN . acadosenv/bin/activate

# # MAKE FILES
# RUN cd acados &&  \
#     make shared_library && \
#     make examples_c && \
#     make run_examples_c

# CMakeFiles
RUN cd acados && ls && \
    mkdir -p build && \
    cd build && \
    cmake .. && \
    make install

RUN pip3 install /root/src/acados/interfaces/acados_template
RUN pip3 install -r /root/src/requirements.txt

CMD ["bash"]

