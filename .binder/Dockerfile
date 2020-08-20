FROM rocker/binder:3.6.3

ARG NB_USER
ARG NB_UID

USER root
COPY . ${HOME}
RUN pwd
RUN ls -la
RUN ls -la .binder 
#RUN R -e 'source("~/.binder/installPackages.R"); installPackages("~/.binder/r-requirements.dcf")'
#RUN R -e 'devtools::install_github("HepingWang/zeptosensPkg", subdir="zeptosensPkg")'
RUN R -e 'devtools::install_deps("zeptosensPkg/zeptosensPkg")'
RUN R -e 'devtools::install("zeptosensPkg/zeptosensPkg")'
RUN chown -R ${NB_USER} ${HOME}

USER ${NB_USER}
