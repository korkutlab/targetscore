# Worked once with CRAN rsconnect_0.8.16.tgz 
install.packages("rsconnect")
#devtools::install_github("rstudio/rsconnect")

shinyapps_user <- Sys.getenv("SHINY_USER")
shinyapps_token <- Sys.getenv("SHINY_TOKEN")
shinyapps_secret <- Sys.getenv("SHINY_SECRET")
rsconnect::setAccountInfo(name=shinyapps_user, token=shinyapps_token, secret=shinyapps_secret)

devtools::install_github('korkutlab/zeptosensPkg', subdir="zeptosensPkg", upgrade='never', dependencies=TRUE)
devtools::install_github('cmap/morpheus.R', upgrade="never")

pkg_dir <- find.package("zeptosensPkg")
shiny_dir <- file.path(pkg_dir, "shiny")

rsconnect::deployApp(appDir=shiny_dir, 
                     appName="targetscore", 
                     appTitle="TargetScore", 
                     forceUpdate=TRUE)
