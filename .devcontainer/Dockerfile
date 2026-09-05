ARG VARIANT=4.5.1
FROM rocker/r-ver:${VARIANT}

RUN useradd -s /bin/bash -m rstudio \
    && usermod -a -G staff rstudio
