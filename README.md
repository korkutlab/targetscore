# TargetScore

TargetScore is a data-driven network modeling algorithm for drug response analysis in cancer.

## Quick Start

Install [Git](https://git-scm.com/), R 3.5 or later, and Pandoc. Pandoc is
included with RStudio; the Docker validation uses R 4.6.

1. Clone the repository and enter its root directory:

```sh
git clone https://github.com/korkutlab/targetscore.git
cd targetscore
```

2. From the repository root, start R and install the package dependencies and
   TargetScore:

```r
install.packages(
  c("remotes", "rmarkdown"),
  repos = "https://cloud.r-project.org"
)
remotes::install_github(
  "cmap/morpheus.R",
  dependencies = TRUE,
  upgrade = "never"
)
remotes::install_deps(
  "targetscore",
  dependencies = TRUE,
  upgrade = "never"
)
remotes::install_local(
  "targetscore",
  dependencies = FALSE,
  upgrade = "never"
)
```

3. In the same R session, render the tutorial:

```r
dir.create("output", showWarnings = FALSE)
output_dir <- normalizePath("output")

rmarkdown::render(
  "targetscore/vignettes/target_score_tutorial.Rmd",
  output_format = "html_document",
  output_dir = output_dir,
  params = list(output_dir = output_dir)
)
```

4. Open `output/target_score_tutorial.html` in a web browser. All other files
   created by the tutorial are also retained in `output/`.

Installation requires access to CRAN and the public `cmap/morpheus.R` GitHub
repository. Some R dependencies may require platform-specific compilers or
system libraries when a prebuilt binary is unavailable.

## Docker Validation and Server

A normal build runs both the Shiny startup check and the complete vignette pipeline before creating the production image:

```sh
docker build --tag targetscore .
docker run --rm --publish 3838:3838 targetscore
```

Open <http://127.0.0.1:3838/> to use the application. The optional `PORT` environment variable changes the container's listening port; publish the same port when overriding it.

Build only through the verification stage:

```sh
docker build --target verify --tag targetscore-verify .
```

Force Buildx to rerun that stage without invalidating the dependency-installation cache:

```sh
docker buildx build --target verify --no-cache-filter verify --load --tag targetscore-verify .
```

Building requires network access to CRAN and the public `cmap/morpheus.R` GitHub repository. Running the production image requires no credentials, mounted data, external services, or network access.
