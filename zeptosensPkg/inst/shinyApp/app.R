#####################################################
#######       Update Module       ###################
#####################################################

# library(rsconnect)
# rsconnect::setAccountInfo(name='hithut', token='B512E4272720432CD9C6F5E47388416D', secret='/Y5XUfYMz02GGocWAcE/vQ+AsSrhpfEA8mBSk8Gl')
# rsconnect::deployApp()
# #

#####################################################
########           UI            ####################
#####################################################
library(pheatmap)
library(shiny)
library(glasso)
library(zeptosensPkg)
library(zeptosensUtils)
library(ggplot2)
library(ggrepel)

ui <-navbarPage(
  #Theme Add In
  #theme="bootstrap.min.css",
  #Or
  #tags$head(
  #tags$link(rel="stylesheet",type="text/css",href="file.css")
  #)
  #Or
  #tags$style(HTML("p{color:red;}"))
  
  #Header
  "Target Score",
  tabPanel("Get Start",
           sidebarPanel(img(width=450,src="Intro.jpg"),
                        tags$small(
                          "Source: TargetScore: ",
                          " July 10, 2005 show at the ",
                          "Heping Wang, and others",
                          a(href="http://commons.wikimedia.org/wiki/User:Sfoskett",
                            "Ref:TargetScore")
                          )),
           mainPanel(
             tags$h2("Why we developed TargetScore"),
             #tags$p("Background",style="font-family:Impact"),
             p("Targeted therapies have been substantially successful in treatment of diverse cancer types
               (Gu et al., 2016; Manzano et al., 2016; Mayekar and Bivona, 2017; Ohmoto and Yachida, 2017;
               Tapia Rico et al., 2017).However, resistance to therapy is virtually inevitable and can manifest 
               as a lack of response to therapy (intrinsic)or disease progression after temporary response 
               (acquired resistance) (Holohan et al., 2013).A recurrent mechanism of resistance is activation
               of compensatory oncogenic pathways (e.g., via feedback loops in short-term or secondary oncogenic
               alterations in long-term) in response to targeting a genomic aberration (Holohan et al., 2013;
               Niederst and Engelman, 2013)A relatively simple way to interrogate adaptive responses to targeted 
               agents is to rank changes in mRNA and (phospho)protein expression for individual genes based on 
               high throughput omics data.  However, the rank-based approach cannot capture collective changes
               that lead to robust phenotypic transitions such as drug resistance.It is also more likely to 
               detect druggable targets within collectively functioning network modules compared to individual 
               genes/proteins. Here, with this motivation, we developed a method to analyze drug response, 
               collective pathway adaptation mechanisms and discover effective drug combinations. Shown as in Fig."),
             
             tags$h2("What does TargetScore do?"),
             #tags$p("Introduction",style="font-family:Impact"),
             p("The algorithm determines the network level, collective responses to perturbations  with a multi-step approach. 
               Target Score Basically have two steps."),
             p("Step1. Contructing Reference Network",style="font-family:Impact"),
             p("In the first step, we generate a reference network model that capture the
               signaling interactions underlying the drug responses. The reference network is inferred using either a knowledge
               based approach with pathway information from multiple signaling databases or a data-driven approach that uses 
               Hybridh Adjusted-Glasso Algorithm (HAGA).There are three provided ways of Contructing Reference Network."),
             strong("predictBionetwork Module:",style = "font-family: 'times'; font-si16pt"),
             p("Contruct the network through Pathway Commons database.Of interest to this study, biochemical reactions, 
               posttranslational modifications, and complex formation, which are all encoded in BioPAX language are taken 
               into account."),
             a(href="https://www.pathwaycommons.org","Pathway Commons"),
             p("The pathway interactions that match to the proteomic species profiled by the RPPA data are
               defined with 4 relation types: phosphorylation, dephosphorylation, expression upregulation, and expression 
               downregulation. The list of gene names and corresponding posttranslational modifications (phosphorylation) 
               is provided to the TargetScore package. The interactions between the molecules of interest are extracted 
               using the signedPC module in BioPAX.The extracted interaction set is evaluated with expert curation and 
               serve as the reference network."),
             strong("predictDatnetwork Module",style = "font-family: 'times'; font-si16pt"),
             p("Data driven reference network captures the molecular associations and inter-tumor heterogeneity across 
               a population of samples with shared characteristics (e.g., ovarian cancer patients).Proteomic datasets, 
               which capture the signaling co-variations serve as  the experimental constraint for network inference. 
               Such datasets can be publicly available (e.g., TCGA data) or custom generated (drug response data, Korkut et al, Elife).
               We inferred from glasso algorithm and provided here as the data driven algorithm in constructing network."),
             a(href="http://statweb.stanford.edu/~tibs/ftp/glasso-bio.pdf","Graphical Lasso algorithm"),
             p("glasso generates a partial correlation based network by estimating a sparse inverse of the covariance matrix 
               using an L1 penalty.The algorithm estimates the precision (i.e., inverse covariance) matrix 
               through maximization of the log-likelihood function by applying the L1 penalty"),
             strong("predictHybnetwork",style = "font-family: 'times'; font-si16pt"),
             p("We have modified the glasso algorithm and developed the Adjusted-glasso algorithm to infer the reference
               network model using priors and with directionality information. Adjusted-glasso algorithm introduces the 
               biology prior information from Pathway Commons Database.The introduction of this symmetric prior information
               matrix can be seen as an adjustment to the inference penalty and allows us to apply varying amounts of 
               penalties to different elements in precision matrix based on prior interactions.The developed algoritm
               give larger probability for the pre-existing/experimentally validated interactions from Pathway Commons database while
               take disease specific global signaling into account."),
             p("Step2. Calculating Target Score Through Reference Network",style="font-family:Impact"),
             p("target score (TS) that quantifies the adaptive pathway responses to a perturbation as a sum of the response from 
               each individual phosphoprotein level and its pathway neighborhood is calculated for each protein in each sample.
               The calculation combines the cell type-specific drug response data with the reference network model information. 
               High target score identifies genes involved in adaptive response (e.g., upregulation of RTK mRNA expression by 
               MEK inhibitor via a feedback loop [[CITE]]) and low target score corresponds to the immediate impact of the drug.")
             
             #tags$a(href="http://www.git.com","TargetScore Pakage Bio"),
           )),
  navbarMenu("Data",
             tabPanel("AntibodyMapFile",
                      wellPanel(
                        h1("AntibodyMapFile"),
                        h2("File Descriptions"),
                        p("AntibodyMapFile is the Gene Name list Database developed by "),
                        a(href="https://odin.mdacc.tmc.edu/~akorkut/#/home","Anil Korkut Group"),
                        p("which provided the onsite AnitobodyLabel of local upload files with the Gene_Symbol used within the Database.While
                          AntibodyMapFile also provided the information of phosphorylation activation/deactivation."),
                        p(""),
                        strong("AntibodyLabel"),
                        p("The AntibodyLabel is the label of local data for each Antibody. The AnitibodyLabel will serve as the label alongside all calculation."),
                        strong("Source"),
                        p("Source of the Antibody. The provided AntibodyMapFile contains two sources including MDACC as MD Anderson Cancer Center
                          and MSKCC standing for Memorial Sloan Kettering Cancer Center."),
                        strong("NodeName"),
                        p("MDACC standardized Antibody Name."),
                        strong("Gene_Symbol"),
                        p("Corresponding HGNC symbol and Ensembl ID"),
                        strong("Sites"),
                        p("phosphorylation site. na stands for no phosphorylation. One site phosphorylation, for example:S473 for Akt_pS473,
                          Two site phosphorylation, for example :Y1234|Y1235 for c.Met_pY1234_Y1235"),
                        strong("Effect"),
                        p("The effect of Phosphorylation. Including:"),
                        p("c : no phosphorylation"),
                        p("a : activation"),
                        p("i : inhibition"),
                        h2("Download"),
                        tags$a(href='data/antibodyMapfile.txt', target='blank', 'antibodyMapfile_localfile', download = 'antibodyMapfile.txt'),
                        h2("Data Usage"),
                        strong("Please cite the following when using these data"),
                        p("Anil K.et al ...")
                      )
             ),
             tabPanel("Global Signaling Data",
                      wellPanel(
                      h1("Global Signaling Data"),
                      h2("File Description"),
                      p("Proteomic datasets, which capture the signaling co-variations serve as  the experimental constraint for network
                        inference. Such datasets can be publicly available (e.g., TCGA data) or custom generated 
                        (drug response data, Korkut et al, Elife).The provided example datas"),
                      p("Coloumns: HGNC symbol and Ensembl ID"),
                      p("Rows: Patients Samples"),
                      h2("DownLoad"),
                      p("Here provided an example from Public Database TCGA of Breast Cancer.Protein level expression data for all genes,
                        Log2 transformed."),
                      tags$a(href='data/TCGA-BRCA-L4_1.csv', target='blank', 'TCGA_BRCA_localfile', download = 'TCGA-BRCA-L4_1.csv'),
                      h2("Data Usage"),
                      strong("Please cite the following when using these data"),
                      p("Anil K.et al ...")
             )),
             tabPanel("Functional Score file",
                      wellPanel(
                      h1("Functional Score file"),
                      h2("File Description"),
                      p("A functional role is assigned as a numeric score to proteomic entities. 
                      there is evidence for an entity’s function as an oncogene or tumor suppressor in cancer. 
                      A central basis for this cancer role comes the curated and widely-used COSMIC database.
                      Two different functional score value were assigned: oncogene as +1,and tumor suppressor as -1.
                      it is also possible for users to manually alter these scores referring from literature
                      or through expert editing, if necessary."),
                      p("gene : HGNC symbol and Ensembl ID"),
                      p("fs : Corresponding functional score"),
                      h2("DownLoad"),
                      tags$a(href='data/Cosmic.txt', target='blank', 'FsFile_localfile', download = 'Cosmic.txt'),
                      h2("Data Usage")
              )),
             tabPanel("","")),
  tabPanel("App",
  tags$img(height=100,width=1000,src="https://3c1703fe8d.site.internapcdn.net/newman/gfx/news/hires/2018/nistbuildsst.jpg"),
  #HTML() use html raw code
  #Side Bar for algorithmn choose
  fluidRow(
    #AntibodyMapFile
    column(3,
           wellPanel(
             fileInput("Antibody","Attach your AntibodyMapFile Signaling File (.csv or .txt)",
                       buttonLabel = "Browse...",
                       placeholder = "No file selected",
                       accept = c(
                         "text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
             
             tags$hr(),
             checkboxInput("header1", "Header", TRUE)
           )),
   #Global signaling data
   column(3,
   wellPanel(
    fileInput("SigData","Attach your Global Signaling File (.csv)",
              buttonLabel = "Browse...",
              placeholder = "No file selected",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")),
    
    tags$hr(),
    checkboxInput("header2", "Header", TRUE)
  )),
    #fs functional score
  column(3,
  wellPanel(
    fileInput("FsFile","Attach your Functional Score File (.csv)",
              buttonLabel = "Browse...",
              placeholder = "No file selected",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")),
    tags$hr(),
    checkboxInput("header3", "Header", TRUE)
  )),
    #Drug perturbation data
  column(3,
  wellPanel(
    fileInput("DrugData","Attach your Drug Perturbation Data(.csv)",
              buttonLabel = "Browse...",
              placeholder = "No file selected",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")),
    tags$hr(),
    checkboxInput("header4", "Header", TRUE)
  ))),
  fluidRow(
    #network contruction algorithmn choice: Bio,Dat,Hyb
  wellPanel(
    selectInput("NetworkAlgorithmn",label="Network Construction Algorithmns:",
                choices = c("Hybrid-driven Network construction Algorithmns" = "Hyb",
                            "Biology-inferral Network construction Algorithmns" = "Bio",
                            "Data-driven Network construction Algorithmns" = "Dat"),
                selected = NULL),
  
    #TS calc type choice:
  
    selectInput("TSCalcType","Target Score Calculation type",
                choices = c("LinebyLine" = "LinebyLine",
                            "Pooled" = "Pooled"),
                selected = NULL),
    textInput("filename", "Filename", "file1"),
    
    #volcano plot line choice
    numericInput("Line","Line Number","1"),
    #max distance of protein network 
    numericInput("max_Dist","Maximum Protein Distance","1"),
    
    actionButton("Submit",label = "Submit/Update",icon = NULL,width = NULL)

    )),
  
# Results showing
  mainPanel(
    tags$h2("Guideline:"),
    p("..."),
    tags$a(href="http://www.git.com","TargetScore Pakage Bio"),
    
    #Results showing in Tabs (can use navlistPanel to show on left)
    tabsetPanel(
    tabPanel("Test module",
    dataTableOutput("Test")),
    tabPanel("Functional Score Value",
             dataTableOutput("FS_value")),
    tabPanel("Heatmap",
    #heatmap of the Data
    plotOutput("heatmap",height = 200,width = 1000)),
    
    tabPanel("TargetScore",
             #heatmap of the Data
             plotOutput("TSheat",height = 200,width = 1000)),
    
    tabPanel("VolcanoPlot of TargetScore",
    #Plots
    plotOutput("volcanoplot"))
  ))),
tabPanel("CART Project",
         sidebarPanel(
           sidebarPanel(
             selectInput("CancerType",label="Cancer Type (Disease type):",
                         choices = c("Breast Cancer" = "BRCA",
                                     "Ovarian Cancer"= "OV",
                                     "Melanoma" = "SKCM",
                                     "Prostate carcinoma" = "PRAD",
                                     "Sarcoma" = "SARC",
                                     "Esophageal carcinoma" = "ESCA",
                                     "B-cell Lymphoma" = "DLBC",
                                     "Uterine Endometrial Carcinoma" = "UCEC",
                                     "Uterine Carcinosarcoma" = "UCS",
                                     "Sarcoma" = "SARC",
                                     "Lung adenocarcinoma" = "LUAD",
                                     "rhabdomyosarcoma" = "rhabdomyosarcoma",
                                     "Pancreatic adenocarcinoma" = "PAAD",
                                     "Prostate adenocarcinoma" = "PRAD",
                                     "Glioma" = "GBMLGG",
                                     "Head and Neck squamous cell carcinoma" = "HNSC",
                                     "Colon adenocarcinoma" = "COAD",
                                     "Kidney cell carcinoma" = "KIRP",
                                     "Unknown" = "NA"
                                     ),
                         selected = NULL),
             selectInput("Cellline",label="Cell Line :",
                         choices = c("OV2008" = "OV2008"
                                     
         
                                     ),
                         selected = NULL),
             selectInput("DrugType",label="Drug:",
                         choices = c("MEKi" = "MEK inhibitor",
                                     "PI3Ki"= "PI3K inhibitor",
                                     "AKTi" = "AKT inhibitor"),
                         selected = NULL),
             selectInput("Dependence",label="TIME/ DOSE Dependence",
                         choices = c("Time" = "TIME",
                                     "Dose"= "Dose"),
                         selected = NULL)
           )
         ),
         mainPanel(
           plotOutput("CARTvolcano"),
           plotOutput("CARTheat")
         )),
tabPanel("About",
         mainPanel(
           p("Target Score is an under-development algorithm with current website version of V1.0.")
         ))
)

#####################################################
########           SERVER        ####################
#####################################################

server <-function(input,output, session, stringsAsFactors){
  
  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.
  
  #read in Inputfiles
  
  #built reactive object through update/submit new dataset
  #SigFile <- eventReactive(input$Submit,{input$SigData}) 
  #FsFile <- eventReactive(input$Submit,{input$FsFile}) 
  #DrugFile <- eventReactive(input$Submit,{input$DrugData}) 
  
  #Algorithm <- eventReactive(input$Submit,{input$NetworkAlgorithmn}) 
  #CalcType <- eventReactive(input$Submit,{input$TSCalcType}) 
  
  #load in system file and Default value and the Default File
  
  #file from input
  
  #Drug File
  DrugDat <- reactive({
    DrugFile<-input$DrugData
    if (is.null(DrugFile))
    return(NULL)
    if(input$header4==T){DrugDat=read.csv(DrugFile$datapath, row.names =1 )}
    if(input$header4==F){DrugDat=read.csv(DrugFile$datapath)}
    return(DrugDat)
  })
  
  #AntibodyMap File (Default at System file) 
  AntiDat <- reactive({
    AntibodyMapFile<-input$Antibody
  if (is.null(AntibodyMapFile))
  AntiDat<-system.file("targetScoreData", "antibodyMapfile_08092019.txt", package = "zeptosensPkg")
  if (!is.null(AntibodyMapFile))
  AntiDat <-read.csv(AntibodyMapFile$datapath, header = input$header1,stringsAsFactors = F )
   return(AntiDat)
  })
  
  #FS File
  FSDat <- reactive({
    FSFile<-input$FsFile
  if (is.null(FSFile))
    return(NULL)
    if(input$header3==T){FSDat <-read.csv(FSFile$datapath,row.names = 1,stringsAsFactors = F )}
    if(input$header3==F){FSDat <-read.csv(FSFile$datapath,stringsAsFactors = F )}
  return(FSDat)
  })
  
  #TSCalcType
  TSCalcType <- reactive({
  TSCalcType<-input$TSCalcType
  return(TSCalcType)
  })

  #filename
  fileName<-reactive({
    fileName<-input$filename
    return(fileName)
    })
  
  #network algorithm
  NetworkAlgorithmn <-reactive({
    NetworkAlgorithmn<-input$NetworkAlgorithmn
    return(NetworkAlgorithmn)})
  
  #Global Signaling file
  SigDat <- reactive({
    SigFile<-input$SigData
  if (is.null(SigFile))
    return(NULL)
  SigDat <-read.csv(SigFile$datapath, header = input$header2,stringsAsFactors = F )
  return(SigDat)
  })
  
  #nProt
  nProt<-reactive({
    data<-DrugDat()
    nProt<-ncol(data)
    return(nProt)
    })
  
  #nCond
  nCond<-reactive({
    data<-DrugDat()
    nCond<-nrow(data)
    return(nCond)
  })
  #Line
  nline<-reactive({
    nline<-input$Line
    return(nline)
  })
  
  #maxDist
  maxDist<-reactive({
    maxDist<-input$max_Dist
    return(maxDist)
  })
  
  #choosing the way to construct reference network
NetworkInferred<-reactive({
    NetworkAlgo<-NetworkAlgorithmn()
    DrugData<-DrugDat()
    AntiData<-AntiDat()
    SigData<-SigDat()
    nPro<-nProt()
    maxiDist<-maxDist()
    
  if (NetworkAlgo=="Bio"){
    # reference network
    network=zeptosensPkg:::predictBioNetwork(nProt =nPro,proteomicResponses = DrugData,antibodyMapFile = AntiData,maxDist = maxiDist)
    wk=network$wk
    wks <- network$wks
    dist_ind <- network$dist_ind
    inter <- network$inter
  }
  
  if(NetworkAlgo=="Dat"){
    network=zeptosensPkg:::predictDatNetwork(data =SigData,nProt=nPro,proteomicResponses=DrugData,maxDist = maxiDist)
    wk=network$wk
    wks <- network$wks
    dist_ind <- network$dist_ind
    inter <- network$inter
  }
  
  if(NetworkAlgo=="Hyb"){
    # prior 
    wk=zeptosensPkg:::predictBioNetwork(nProt =nPro,proteomicResponses = DrugData,maxDist = maxiDist,antibodyMapFile = AntiData)
    #Hyb
    network=zeptosensPkg:::predictHybNetwork(data =SigData,prior=wk,nProt=nPro,proteomicResponses=DrugData)
    
    wk=network$wk
    wks <- network$wks
    dist_ind <- network$dist_ind
    inter <- network$inter
  }
    NetworkInferred<-list(wk=wk,wks=wks,dist_ind=dist_ind,inter=inter)
    return(NetworkInferred)
  })
    #fs
FsValue<-reactive({
  DrugData<-DrugDat()
  AntiData<-AntiDat()
  nPro<-nProt()
    if(is.null(input$FsFile))
      {FsValue=zeptosensPkg:::getFsVals(nProt =nPro ,proteomicResponses =DrugData,antibodyMapFile = AntiData)}
    if(!is.null(input$FsFile))
      {FsValue=zeptosensPkg:::getFsVals(nProt =nPro ,proteomicResponses =DrugData,fsValueFile=FSDat,antibodyMapFile = AntiData)}
return(FsValue)
  })

    #Calc TargetScore
TS.r<-reactive({
  
  #call up reactive items
  NetworkInferred<-NetworkInferred()
  DrugDat<-DrugDat()
  nProt<-nProt()
  FsValue<-FsValue()
  fileName<-fileName()
  maxDist<-maxDist()
  
  #give output filename
  targetScoreOutputFile <-paste0(fileName,"TS.txt")
  matrixWkOutputFile <-paste0(fileName,"wk.txt")
  signedMatrixWkOutputFile <-paste0(fileName,"wks.txt")
  
  #Network inferred
  wk=NetworkInferred$wk
  wks <- NetworkInferred$wks
  dist_ind <- NetworkInferred$dist_ind
  inter <- NetworkInferred$inter
  
  if(TSCalcType=="LinebyLine"){
      
      #Calc Std
      DrugDat[is.na(DrugDat)]<-0
      stdev <- zeptosensPkg:::sampSdev(nSample=nrow(DrugDat),nProt=ncol(DrugDat),nDose=1,nX=DrugDat)
      #normalization
      proteomicResponses<- DrugDat
      for(i in 1:nProt){
        for (j in 1:nrow(proteomicResponses)){
          proteomicResponses[j,i] <- (DrugDat[j,i]/stdev[i])      
        }
      }
      
      #Bootstrap in Getting TargetScore
      TS <- array(0,dim=c(nCond,nProt))
      TS.p <- array(0,dim=c(nCond,nProt))
      TS.q <- array(0,dim=c(nCond,nProt))
      
      nPerm=1000
      
      for(i in 1:nCond){
        
        results <- zeptosensPkg:::getTargetScore(wk=wk,
                                                 wks=wks,
                                                 dist_ind=dist_ind,
                                                 inter=inter,
                                                 nDose=1, 
                                                 nProt=nProt, 
                                                 proteomicResponses=proteomicResponses[i,], 
                                                 maxDist=maxDist, 
                                                 nPerm=nPerm,
                                                 cellLine=fileName, 
                                                 verbose=FALSE,fsFile=FsValue,
                                                 targetScoreOutputFile=targetScoreOutputFile, 
                                                 matrixWkOutputFile=matrixWkOutputFile,
                                                 targetScoreQValueFile="q.txt", 
                                                 targetScoreDoseFile="TS_d.txt",
                                                 targetScorePValueFile="p.txt")
        TS[i,]=results$ts
        TS.p[i,]=results$pts
        TS.q[i,]=results$q
      }
      colnames(TS)=colnames(DrugDat)
      rownames(TS)=rownames(DrugDat)

      colnames(TS.p)=colnames(DrugDat)
      rownames(TS.p)=rownames(DrugDat)

      colnames(TS.q)=colnames(DrugDat)
      rownames(TS.q)=rownames(DrugDat) 
  }
  
  if(TSCalcType=="Pooled"){
    
    #Calc Std
    DrugDat[is.na(DrugDat)]<-0
    stdev <- zeptosensPkg:::sampSdev(nSample=nrow(DrugDat),nProt=ncol(DrugDat),nDose=1,nX=DrugDat)
    #normalization
    proteomicResponses<- DrugDat
    for(i in 1:nProt){
      for (j in 1:nrow(proteomicResponses)){
        proteomicResponses[j,i] <- (DrugDat[j,i]/stdev[i])      
      }
    }
    
    nPerm=1000
      
      results <- zeptosensPkg:::getTargetScore(wk=wk,
                                               wks=wks,
                                               dist_ind=dist_ind,
                                               inter=inter,
                                               nDose=1, 
                                               nProt=nProt, 
                                               proteomicResponses=proteomicResponses, 
                                               maxDist=maxDist, 
                                               nPerm=nPerm,
                                               cellLine=fileName, 
                                               verbose=FALSE,fsFile=FsValue,
                                               targetScoreOutputFile=targetScoreOutputFile, 
                                               matrixWkOutputFile=matrixWkOutputFile,
                                               targetScoreQValueFile="q.txt", 
                                               targetScoreDoseFile="TS_d.txt",
                                               targetScorePValueFile="p.txt")
      TS=results$ts
      TS.p=results$pts
      TS.q=results$q
    
    colnames(TS)=colnames(DrugDat)
    
    colnames(TS.p)=colnames(DrugDat)
    
    colnames(TS.q)=colnames(DrugDat)
  }
    TS.r<-list(TS=TS,TS.p=TS.p,TS.q=TS.q)
    return(TS.r)
})

#Data heatmap
output$heatmap <- renderPlot({
  DrugFile <- DrugDat()
  
  maxDat=max(as.matrix(DrugFile))
  minDat=min(as.matrix(DrugFile))
  bk <- c(seq(minDat,-0.01,by=0.01),seq(0,maxDat,by=0.01))
  data=as.matrix(DrugFile)
  pheatmap(data,
           scale = "none",
           color = c(colorRampPalette(colors = c("navy","white"))(length(seq(minDat,-0.01,by=0.01))),colorRampPalette(colors = c("white","firebrick3"))(length(seq(0,maxDat,by=0.01)))),
           legend_breaks=seq(minDat,maxDat,2),cellwidth = 2, cellheight = 2, fontsize=2, fontsize_row=2,
           breaks=bk)
}) 

output$FS_value <- renderDataTable({
  Data <- FSDat()
}) 
###########################################################################################################################
###########                                Test module                                        #############################
###########################################################################################################################
output$Test <- renderDataTable({
  Data <- AntiDat()
}) 


###########################################################################################################################
###########                                network visualization module                       #############################
###########################################################################################################################




###########################################################################################################################
###########                                TS heatmap module                                 #############################
###########################################################################################################################

output$TSheat <- renderPlot({
       TS.r=TS.r()
       TS=TS.r$TS
       maxDat=max(as.matrix(TS))
       minDat=min(as.matrix(TS))
       bk <- c(seq(minDat,-0.01,by=0.01),seq(0,maxDat,by=0.01))
       data=as.matrix(TS)
       pheatmap(data,
                scale = "none",
                color = c(colorRampPalette(colors = c("navy","white"))(length(seq(minDat,-0.01,by=0.01))),colorRampPalette(colors = c("white","firebrick3"))(length(seq(0,maxDat,by=0.01)))),
                legend_breaks=seq(minDat,maxDat,2),cellwidth = 2, cellheight = 2, fontsize=2, fontsize_row=2,
                breaks=bk)
})

###########################################################################################################################
###########                                TS Volcano plot module                             #############################
###########################################################################################################################

output$volcanoplot <- renderPlot({  
  TS.r=TS.r()
  nline=nline()
   TS=TS.r$TS[nline,]
   TS.q=TS$TS.q[nline,]
   TS<- as.matrix(TS)
   Padj<- as.matrix(TS.q)
# 
   if(nrow(Padj)!=nrow(TS)){
     stop("ERROR:Tag of TS and Qvalue does not match.")
   }
   tmpDat <- data.frame(cbind(TS,-1*log10(Padj)))
   colnames(tmpDat) <- c("TS","neglogQ")
#   
   color <- ifelse(Padj>0.4,"not significant","significant")
   rownames(color) <- rownames(TS)
   tmpDat$labelnames <-  row.names(tmpDat)
   sig01 <- subset(tmpDat, tmpDat$neglogQ > -1*log10(0.4))
   siglabel <- sig01$labelnames
   tmpDat$color <- color
#   
   ggplot() +
     geom_point(data=tmpDat, aes(x=TS, y=neglogQ, color=color), alpha=0.4, size=2) +
     theme_bw() +
     xlab("<TS>") + ylab("-log10 (Q-Value)") + ggtitle("")+
     scale_color_manual(name="", values=c("black", "red"))+
     geom_label_repel(data=sig01, aes(x=sig01$TS, y=sig01$neglogQ,label=siglabel), size=5)
 })
}

shinyApp(server = server,ui=ui)

