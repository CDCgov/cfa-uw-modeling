FROM docker.io/rocker/r-ver:4.4.0

RUN echo 'options(repos = c(P3M = "https://packagemanager.posit.co/cran/__linux__/jammy/latest", \
                            CRAN = "https://cloud.r-project.org"))' >>"${R_HOME}/etc/Rprofile.site"

RUN install2.r remotes
RUN installGithub.r statnet/tergmLite statnet/EpiModelHPC
COPY epimodelcfa_*.tar.gz .
RUN R CMD INSTALL epimodelcfa_*.tar.gz --dependencies=TRUE
