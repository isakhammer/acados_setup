# Essential
FROM ubuntu:20.04
RUN apt-get update
RUN apt-get -y upgrade
RUN apt-get install -yq apt-utils dialog
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get install -yq build-essential software-properties-common
RUN apt-get update

# Install bare bons python
RUN apt-get -y install python3
RUN apt-get -y install python3-pip
RUN apt-get install -yq python3-tk

# Deep bug somwhere when installing matplotlib on ubuntu20
# https://stackoverflow.com/questions/25674612/ubuntu-14-04-pip-cannot-upgrade-matplotllib
RUN apt-get install -yq libfreetype6-dev libxft-dev

RUN /usr/bin/python3 -m pip install --upgrade pi
#RUN pip3 install --upgrade pip

# Install python dependencies
RUN mkdir -p /root/src/
WORKDIR /root/src
COPY . .

# CMake
RUN apt-get install cmake -y


#Acados
RUN cd acados && ls && \
    mkdir -p build && \
    cd build && \
    cmake -DACADOS_WITH_QPOASES=ON ACADOS_PYTHON=ON ACADOS_EXAMPLES=ON \
            -DACADOS_INSTALL_DIR=/root/src/acados/ .. && \
    make install  && \
    cd .. && \
    pip3 install interfaces/acados_template

RUN echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/root/src/acados/lib"" >> /root/.bashrc
RUN echo "export ACADOS_SOURCE_DIR="/root/src/acados"" >> /root/.bashrc


# RUN pip3 install -e /root/src/acados/interfaces/acados_template

# RUN pip3 install -r requirements.txt
