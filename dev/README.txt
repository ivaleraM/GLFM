The files inside these folders are used for development, to test installation procedures inside a docker (See https://docs.docker.com/ for further documentation)

Steps to be able to test with docker:

1- Install docker (stable version, see url above)

2- Run commands:

        sudo docker build -t custom/glfm -f Dockerfile_python .
        sudo docker run -ti --rm custom/glfm

(If you are testing installation for R, write Dockerfile_R instead of
Dockerfile_python)

An independent docker environment will be created with the a priori
requirements already available. The only remaining thing to do is to test the install procedure:

For PYTHON:

- Clone GLFM git repository:
    git clone https://github.com/ivaleraM/GLFM.git

- Go to install/ and execute:
    bash install_for_python.sh
