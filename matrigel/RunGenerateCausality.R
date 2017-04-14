library(rmarkdown)

outputDir <- "./"
dataDir <- "./hcc1954"

files <- dir(path=dataDir, pattern="^TS_.*.txt$")
dataFileBases <- gsub("\\.txt", "", files)
runHead <- c(TRUE, FALSE)   
datRows <- 10
geneCentric <- c(TRUE, FALSE)

opts <- expand.grid(dataFileBase=dataFileBases, runHead=runHead, datRows=datRows, geneCentric=geneCentric, stringsAsFactors=FALSE)

errors <- NULL

for(i in 1:nrow(opts)) {
  cat("I: ", i, "\n")
  tmpOpts <- opts[i,]
  
  tryCatch({
    rmarkdown::render(file.path("causality.Rmd"), 
                      params = list(dataFileBase=tmpOpts$dataFileBase, runHead=tmpOpts$runHead, 
                                    datRows=tmpOpts$datRows, dataDir=dataDir, geneCentric=tmpOpts$geneCentric), 
                      output_dir = outputDir)    
  }, error=function(e) {
    print(paste('ERROR: ',e))
    
    errors <- rbind(errors, tmpOpts)
  })
}
