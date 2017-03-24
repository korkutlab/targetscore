library(rmarkdown)

outputDir <- "./"

files <- dir(pattern = "^TS_.*.txt$")
dataFileBases <- gsub("\\.txt", "", files)
runHead <- c(TRUE, FALSE)   
datRows <- 10
graphTypes <- c("compatible", "conflicting")

opts <- expand.grid(dataFileBase=dataFileBases, runHead=runHead, datRows=datRows, graphType=graphTypes, stringsAsFactors=FALSE)

errors <- NULL

for(i in 1:nrow(opts)) {
  cat("I: ", i, "\n")
  tmpOpts <- opts[i,]
  
  tryCatch({
    rmarkdown::render(file.path("causalitySphereoids.Rmd"), 
                      params = list(dataFileBase=tmpOpts$dataFileBase, runHead=tmpOpts$runHead, 
                                    datRows=tmpOpts$datRows, graphType=tmpOpts$graphType), 
                      output_dir = outputDir)    
  }, error=function(e) {
    errors <- rbind(errors, tmpOpts)
  })
}
