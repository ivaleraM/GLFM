# Use an official Python runtime as a parent image
FROM ubuntu:xenial

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

# Install Git and other basic modules
RUN apt-get update && \
    apt-get -y --no-install-recommends install \
    ca-certificates \
    git \
    wget \
    bzip2 \
    vim \
    build-essential

# install ANACONDA platform

RUN useradd -m -s /bin/bash anaconda
USER anaconda

# Set the working directory
WORKDIR /home/anaconda

RUN echo 'export PATH=/home/anaconda/conda/bin:$PATH' > /home/anaconda/.bashrc && \
    wget --quiet https://repo.continuum.io/archive/Anaconda2-4.4.0-Linux-x86_64.sh -O /home/anaconda/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /home/anaconda/conda && \
    rm ~/anaconda.sh

# install QT dependencies for plotting functionality
USER root
RUN apt-get -y install python-qt4

# ---------------
# executable part
# ---------------
#
USER anaconda
RUN git clone https://github.com/ivaleraM/GLFM.git
RUN /home/anaconda/GLFM/install/install_for_python.sh

# Run app.py when the container launches
CMD ["python", "/home/anaconda/GLFM/demos/demo_GLFM_infer_simple_call.py"]

#ADD install.sh /home/anaconda
#RUN mv install.sh GLFM/install.sh
#WORKDIR /home/anaconda/GLFM
#RUN bash install.sh

# Copy the current directory contents into the container at /app
#ADD . /app

# Install any needed packages specified in requirements.txt
#RUN pip install -r requirements.txt

# Make port 80 available to the world outside this container
#EXPOSE 80

# Define environment variable
#ENV NAME World

# Run app.py when the container launches
#CMD ["python", "app.py"]
