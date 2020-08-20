FROM rocker/binder:3.6.3

ARG NB_USER
ARG NB_UID

USER root
COPY . ${HOME}
RUN ls -la 
RUN R -e 'source("~/.binder/installPackages.R")'
RUN chown -R ${NB_USER} ${HOME}

USER ${NB_USER}
