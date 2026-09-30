# TargetScore
TargeScore is a data-driven network modeling algorithm for drug response analysis in cancer.

# Quickstart

```
# INSTALL ---- 
install.packages(c("remotes", "rmarkdown")) # remotes to install, rmarkdown to run tutorial RMarkdown file

remotes::install_github("korkutlab/targetscore", subdir="targetscore", upgrade_dependencies=FALSE)
 
# RUN ----
# Download code from GitHub (includes sample data)
download.file(url="https://github.com/korkutlab/targetscore/archive/master.zip", destfile="~/targetscore-master.zip")

# Unzip the .zip file
unzip(zipfile="~/targetscore-master.zip", exdir="~/")

# Switch to the targetscore-master directory
setwd(dir = "~/targetscore-master/targetscore")

# Delete any existing output and make the output directory
output_dir <- "~/targetscore-master/output"
if(dir.exists(output_dir)) {
  unlink(output_dir)
}
dir.create(output_dir)

# Render the RMarkdown tutorial
rmarkdown::render("vignettes/target_score_tutorial.Rmd", "html_document", params = list(output_dir=output_dir))

# The rendered HTML file will be in the same directory as the RMarkdown file 
# Also various outputs will be generated at "~/targetscore-master/output"
```
