ARG R_VERSION=4.6

FROM rocker/r-ver:${R_VERSION} AS base

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        curl \
        git \
        libcurl4-openssl-dev \
        libssl-dev \
        libuv1-dev \
        libxml2-dev \
        pandoc \
        pkg-config \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install dependencies before copying the package source to preserve this cache
# when only project code changes.
RUN Rscript -e 'install.packages(c("data.table", "DT", "ggplot2", "ggrepel", "glasso", "gplots", "knitr", "markdown", "matrixStats", "pheatmap", "plotly", "plyr", "RColorBrewer", "remotes", "rmarkdown", "shiny"), repos = "https://cloud.r-project.org", Ncpus = parallel::detectCores())' \
    && Rscript -e 'remotes::install_github("cmap/morpheus.R", dependencies = TRUE, upgrade = "never")'

WORKDIR /opt/targetscore
COPY targetscore/ ./

RUN R CMD INSTALL .


FROM base AS verify

COPY docker_verify.sh /usr/local/bin/verify-targetscore
RUN /usr/local/bin/verify-targetscore


# This optional target renders a copy of the vignette into a mounted /output
# directory. The build-time verification output remains isolated in `verify`.
FROM verify AS vignette

RUN mkdir -p /output
VOLUME ["/output"]

CMD ["Rscript", "-e", "rmarkdown::render('/opt/targetscore/vignettes/target_score_tutorial.Rmd', output_format = 'html_document', output_dir = '/output', params = list(output_dir = '/output'))"]


FROM verify AS runtime-libraries

RUN rm -rf \
        /usr/local/lib/R/site-library/knitr \
        /usr/local/lib/R/site-library/remotes \
        /usr/local/lib/R/site-library/rmarkdown


FROM rocker/r-ver:${R_VERSION} AS production

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        libcurl4t64 \
        libssl3t64 \
        libuv1t64 \
        libxml2 \
    && rm -rf /var/lib/apt/lists/*

# Copy only the verified installed R libraries. Source, logs, build tools, and
# pipeline artifacts are not included in the production image.
COPY --from=runtime-libraries /usr/local/lib/R/site-library /usr/local/lib/R/site-library

ENV PORT=3838
EXPOSE 3838

CMD ["Rscript", "-e", "shiny::runApp(system.file('shiny', package = 'targetscore'), host = '0.0.0.0', port = as.integer(Sys.getenv('PORT', '3838')))"]
