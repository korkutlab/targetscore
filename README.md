# Installation and Usage 

    install.packages("devtools")
    
    library(devtools)
    install_bitbucket("cbio_mskcc/kirc_immune_project",
        auth_user="discoverUser",
        password="discoverUserPassword",
        build_vignette=TRUE,
        dependencies=TRUE,
        args="--no-multiarch")
        
    library(kircImmuneProject)
    runShinyApp()
