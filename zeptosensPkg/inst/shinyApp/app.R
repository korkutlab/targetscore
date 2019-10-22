#####################################################
#######       Update Module       ###################
#####################################################

# library(rsconnect)
# rsconnect::setAccountInfo(name='hithut', token='B512E4272720432CD9C6F5E47388416D',
#   secret='/Y5XUfYMz02GGocWAcE/vQ+AsSrhpfEA8mBSk8Gl')
# rsconnect::deployApp()

#####################################################
########           UI            ####################
#####################################################
library(pheatmap)
library(shiny)
library(glasso)
library(zeptosensPkg)
library(ggplot2)
library(ggrepel)
library(markdown)

ui <- navbarPage(
  # Theme Add In
  # theme="bootstrap.min.css",

  # Header
  "Target Score",
  tabPanel(
    "Get Start",
    sidebarPanel(
      img(width = 450, src = "Intro.jpg"),
      tags$small(
        "Source: Targetscore: ",
        " July 10, 2005 show at the ",
        "Heping Wang, and others",
        a(
          href = "http://commons.wikimedia.org/wiki/User:Sfoskett",
          "Ref:Targetscore"
        )
      )
    ),
    mainPanel(
      includeMarkdown("www/ts_intro_p1.md"),
      includeMarkdown("www/ts_intro_p2.md")
    )
  ),
  navbarMenu(
    "data",
    tabPanel(
      "antibody_map_file",
      wellPanel(
        h1("antibody_map_file"),
        h2("File Descriptions"),
        p("antibody_map_file is the Gene Name list database developed by "),
        a(href = "https://odin.mdacc.tmc.edu/~akorkut/#/home", "Anil Korkut Group"),
        p("which provided the onsite AnitobodyLabel of local upload files with the Gene_Symbol used within the database.
        While antibody_map_file also provided the information of phosphorylation activation/deactivation."),
        p(""),
        strong("AntibodyLabel"),
        p("The AntibodyLabel is the label of local data for each Antibody. The AnitibodyLabel will serve as the label 
          alongside all calculation."),
        strong("Source"),
        p("Source of the Antibody. The provided antibody_map_file contains two sources including MDACC as MD Anderson 
          Cancer Center and MSKCC standing for Memorial Sloan Kettering Cancer Center."),
        strong("NodeName"),
        p("MDACC standardized Antibody Name."),
        strong("Gene_Symbol"),
        p("Corresponding HGNC symbol and Ensembl ID"),
        strong("Sites"),
        p("phosphorylation site. na stands for no phosphorylation. One site phosphorylation, for example:S473 for 
        Akt_pS473, Two site phosphorylation, for example :Y1234|Y1235 for c.Met_pY1234_Y1235"),
        strong("Effect"),
        p("The effect of Phosphorylation. Including:"),
        p("c : no phosphorylation"),
        p("a : activation"),
        p("i : inhibition"),
        h2("Download"),
        tags$a(
          href = "data/antibodyMapfile.txt", target = "blank", "antibodyMapfile_localfile",
          download = "antibodyMapfile.txt"
        ),
        h2("data Usage"),
        strong("Please cite the following when using these data"),
        p("Anil K.et al ...")
      )
    ),
    tabPanel(
      "Global signaling data",
      wellPanel(
        h1("Global signaling data"),
        h2("File Description"),
        p("Proteomic datasets, which capture the signaling co-variations serve as  the experimental constraint for 
          network inference. Such datasets can be publicly available (e.g., TCGA data) or custom generated (drug 
          response data, Korkut et al, Elife).The provided example datas"),
        p("Coloumns: HGNC symbol and Ensembl ID"),
        p("Rows: Patients Samples"),
        h2("DownLoad"),
        p("Here provided an example from Public database TCGA of Breast Cancer.Protein level expression data for all 
          genes, Log2 transformed."),
        tags$a(
          href = "data/TCGA-BRCA-L4_1.csv", target = "blank", "TCGA_BRCA_localfile",
          download = "TCGA-BRCA-L4_1.csv"
        ),
        h2("data Usage"),
        strong("Please cite the following when using these data"),
        p("Anil K.et al ...")
      )
    ),
    tabPanel(
      "Functional Score file",
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
        tags$a(href = "data/Cosmic.txt", target = "blank", "fs_file_localfile", download = "Cosmic.txt"),
        h2("data Usage")
      )
    ),
    tabPanel("", "")
  ),
  tabPanel(
    "App",
    tags$img(
      height = 100, width = 1000,
      src = "https://3c1703fe8d.site.internapcdn.net/newman/gfx/news/hires/2018/nistbuildsst.jpg"
    ),
    # HTML() use html raw code
    # Side Bar for algorithmn choose
    fluidRow(
      # antibody_map_file
      column(
        3,
        wellPanel(
          fileInput("Antibody", "Attach your antibody_map_file signaling File (.csv or .txt)",
            buttonLabel = "Browse...",
            placeholder = "No file selected",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"
            )
          ),

          tags$hr(),
          checkboxInput("header1", "Header", TRUE)
        )
      ),
      # Global signaling data
      column(
        3,
        wellPanel(
          fileInput("sig", "Attach your Global signaling File (.csv)",
            buttonLabel = "Browse...",
            placeholder = "No file selected",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"
            )
          ),

          tags$hr(),
          checkboxInput("header2", "Header", TRUE)
        )
      ),
      # fs functional score
      column(
        3,
        wellPanel(
          fileInput("fs_file", "Attach your Functional Score File (.csv)",
            buttonLabel = "Browse...",
            placeholder = "No file selected",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"
            )
          ),
          tags$hr(),
          checkboxInput("header3", "Header", TRUE)
        )
      ),
      # Drug perturbation data
      column(
        3,
        wellPanel(
          fileInput("drug_data", "Attach your Drug Perturbation data(.csv)",
            buttonLabel = "Browse...",
            placeholder = "No file selected",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"
            )
          ),
          tags$hr(),
          checkboxInput("header4", "Header", TRUE)
        )
      )
    ),
    fluidRow(
      # network contruction algorithmn choice: Bio,Dat,Hyb
      wellPanel(
        selectInput("network_algorithm",
          label = "Network Construction Algorithmns:",
          choices = c(
            "Hybrid-driven Network construction Algorithmns" = "Hyb",
            "Biology-inferral Network construction Algorithmns" = "Bio",
            "data-driven Network construction Algorithmns" = "Dat"
          ),
          selected = NULL
        ),

        # ts calc type choice:

        selectInput("tsCalcType",
          label = "Target Score Calculation type:",
          choices = c(
            "Line by Line" = "lnl",
            "Pooled" = "pop"
          ),
          selected = NULL
        ),
        textInput("filename", "Filename", "file1"),

        # volcano plot line choice
        numericInput("Line", "Line Number", "1"),
        # max distance of protein network
        numericInput("max_Dist", "Maximum Protein Distance", "1"),

        actionButton("Submit", label = "Submit/Update", icon = NULL, width = NULL)
      )
    ),

    # Results showing
    mainPanel(
      tags$h2("Guideline:"),
      p("..."),
      tags$a(href = "http://www.git.com", "Targetscore Pakage Bio"),

      # Results showing in Tabs (can use navlistPanel to show on left)
      tabsetPanel(
        tabPanel(
          "test module",
          textOutput("test")
        ),

        tabPanel(
          "Antibody Map",
          dataTableOutput("anti_map")
        ),

        tabPanel(
          "Functional Score Value",
          dataTableOutput("fs_value")
        ),

        tabPanel(
          "Heatmap",
          # heatmap of the data
          plotOutput("heatmap", height = 200, width = 1000)
        ),

        tabPanel(
          "edgelist of Network",
          # edgelist of the data
          dataTableOutput("edgelist")
        ),

        tabPanel(
          "Targetscore",
          # heatmap of the data
          plotOutput("tsheat", height = 200, width = 1000)
        ),

        tabPanel(
          "VolcanoPlot of Targetscore",
          # Plots
          plotOutput("volcanoplot")
        )
      )
    )
  ),
  tabPanel(
    "CART Project",
    sidebarPanel(
      sidebarPanel(
        selectInput("CancerType",
          label = "Cancer Type (Disease type):",
          choices = c(
            "Breast Cancer" = "BRCA",
            "Ovarian Cancer" = "OV",
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
          selected = NULL
        ),
        selectInput("Cellline",
          label = "Cell Line :",
          choices = c("OV2008" = "OV2008"),
          selected = NULL
        ),
        selectInput("DrugType",
          label = "Drug:",
          choices = c(
            "MEKi" = "MEK inhibitor",
            "PI3Ki" = "PI3K inhibitor",
            "AKTi" = "AKT inhibitor"
          ),
          selected = NULL
        ),
        selectInput("Dependence",
          label = "TIME/ DOSE Dependence",
          choices = c(
            "Time" = "TIME",
            "Dose" = "Dose"
          ),
          selected = NULL
        )
      )
    ),
    mainPanel(
      plotOutput("CARTvolcano"),
      plotOutput("CARTheat")
    )
  ),
  tabPanel(
    "About",
    mainPanel(
      p("Target Score is an under-development algorithm with current website version of V1.0.")
    )
  )
)

#####################################################
########           SERVER        ####################
#####################################################

server <- function(input, output, session, strings_as_factors) {

  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.

  # read in Inputfiles

  # built reactive object through update/submit new dataset
  # sig_file <- eventReactive(input$Submit,{input$sig_data})
  # fs_file <- eventReactive(input$Submit,{input$fs_file})
  # drug_file <- eventReactive(input$Submit,{input$drug_data})

  # Algorithm <- eventReactive(input$Submit,{input$network_algorithm})
  # CalcType <- eventReactive(input$Submit,{input$tsCalcType})

  # load in system file and Default value and the Default File

  # file from input

  # Drug File
  drug_dat <- reactive({
    drug_file <- input$drug_data
    if (is.null(drug_file)) {
      return(NULL)
    }
    if (input$header4 == TRUE) {
      drug_dat <- read.csv(drug_file$datapath, row.names = 1)
    }
    if (input$header4 == FALSE) {
      drug_dat <- read.csv(drug_file$datapath)
    }
    return(drug_dat)
  })

  # AntibodyMap File (Default at System file)
  anti_dat <- reactive({
    antibody_map_file <- input$Antibody
    if (is.null(antibody_map_file)) {
      anti_dat <- system.file("targetScoreData", "antibodyMapfile_08092019.txt", package = "zeptosensPkg")
    }
    if (!is.null(antibody_map_file)) {
      anti_dat <- read.csv(antibody_map_file$datapath, header = input$header1, strings_as_factors = FALSE)
    }
    return(anti_dat)
  })

  # FS File
  fs_dat <- reactive({
    fs_file <- input$fs_file
    if (is.null(fs_file)) {
      return(NULL)
    }
    if (input$header3 == TRUE) {
      fs_dat <- read.csv(fs_file$datapath, header = TRUE, strings_as_factors = FALSE)
    }
    if (input$header3 == FALSE) {
      fs_dat <- read.csv(fs_file$datapath, header = FALSE, strings_as_factors = FALSE)
    }
    return(fs_dat)
  })

  # tsCalcType
  ts_type <- reactive({
    ts_type <- input$tsCalcType
    return(ts_type)
  })

  # filename
  file_name <- reactive({
    file_name <- input$filename
    return(file_name)
  })

  # network algorithm
  network_algorithm <- reactive({
    network_algorithm <- input$network_algorithm
    return(network_algorithm)
  })

  # Global signaling file
  sig_dat <- reactive({
    sig_file <- input$sig
    if (is.null(sig_file)) {
      return(NULL)
    }
    if (input$header2 == TRUE) {
      sig_dat <- read.csv(sig_file$datapath, row.names = 1, strings_as_factors = FALSE)
    }
    if (input$header2 == FALSE) {
      sig_dat <- read.csv(sig_file$datapath, header = FALSE, strings_as_factors = FALSE)
    }
    return(sig_dat)
  })

  # n_prot
  n_prot <- reactive({
    data <- drug_dat()
    n_prot <- ncol(data)
    return(n_prot)
  })

  # n_cond
  n_cond <- reactive({
    data <- drug_dat()
    n_cond <- nrow(data)
    return(n_cond)
  })
  # Line
  nline <- reactive({
    nline <- input$Line
    return(nline)
  })

  # max_dist
  max_dist <- reactive({
    max_dist <- input$max_Dist
    return(max_dist)
  })

  # choosing the way to construct reference network
  network_inferred <- reactive({
    network_algo <- network_algorithm()
    drug_data <- drug_dat()
    anti_data <- anti_dat()
    sig_data <- sig_dat()
    n_pro <- n_prot()
    maxi_dist <- max_dist()

    if (network_algo == "Bio") {
      # reference network
      network <- zeptosensPkg::predict_bio_network(
        n_prot = n_pro,
        proteomic_responses = drug_data,
        antibody_map_file = anti_data,
        max_dist = maxi_dist
      )
      wk <- network$wk
      wks <- network$wks
      dist_ind <- network$dist_ind
      inter <- network$inter
    }

    if (network_algo == "Dat") {
      network <- zeptosensPkg::predict_dat_network(
        data = sig_data,
        n_prot = n_pro,
        proteomic_responses = drug_data,
        max_dist = maxi_dist
      )
      wk <- network$wk
      wks <- network$wks
      dist_ind <- network$dist_ind
      inter <- network$inter
    }

    if (network_algo == "Hyb") {
      # prior
      wk <- zeptosensPkg::predict_bio_network(
        n_prot = n_pro,
        proteomic_responses = drug_data,
        max_dist = maxi_dist,
        antibody_map_file = anti_data
      )$wk
      # Hyb
      network <- zeptosensPkg::predict_hyb_network(
        data = sig_data,
        prior = wk,
        n_prot = n_pro,
        proteomic_responses = drug_data
      )

      wk <- network$wk
      wks <- network$wks
      dist_ind <- network$dist_ind
      inter <- network$inter
    }
    network_inferred <- list(wk = wk, wks = wks, dist_ind = dist_ind, inter = inter)
    return(network_inferred)
  })
  # fs
  fs_value <- reactive({
    drug_data <- drug_dat()
    anti_data <- anti_dat()
    n_pro <- n_prot()
    fs_dat <- fs_dat()
    if (is.null(fs_dat)) {
      fs_value <- zeptosensPkg::get_fs_vals(
        n_prot = n_pro,
        proteomic_responses = drug_data,
        antibody_map_file = anti_data
      )
    }
    if (!is.null(fs_dat)) {
      fs_value <- zeptosensPkg::get_fs_vals(
        n_prot = n_pro,
        proteomic_responses = drug_data,
        fs_valueFile = fs_dat,
        antibody_map_file = anti_data
      )
    }
    return(fs_value)
  })

  # Calc Targetscore
  ts_r <- reactive({

    # call up reactive items
    ts_type <- ts_type()
    network_inferred <- network_inferred()
    drug_data <- drug_dat()
    n_pro <- n_prot()
    n_con <- n_cond()
    fs_value <- fs_dat() ### Warning :#Shoule be changed back to fs_value when fs_value() Reactive Back to Work ###
    file_name <- file_name()
    maxi_dist <- max_dist()

    # Network inferred
    wk <- network_inferred$wk
    wks <- network_inferred$wks
    dist_ind <- network_inferred$dist_ind
    inter <- network_inferred$inter

    if (ts_type == "lnl") {
      drug_data[is.na(drug_data)] <- 0

      # Calc Std (Normalization request)
      # stdev <- zeptosensPkg::samp_sdev(nSample=nrow(drug_data),n_prot=ncol(drug_data),n_dose=1,nX=drug_data)
      # #normalization
      # proteomic_responses<- drug_data
      # for(i in 1:n_prot) {
      #   for (j in 1:nrow(proteomic_responses)) {
      #     proteomic_responses[j,i] <- (drug_data[j,i]/stdev[i])
      #   }
      # }
      #
      proteomic_responses <- drug_data
      # Bootstrap in Getting Targetscore
      ts <- array(0, dim = c(n_con, n_pro))
      ts_p <- array(0, dim = c(n_con, n_pro))
      ts_q <- array(0, dim = c(n_con, n_pro))

      n_perm <- 1000

      for (i in 1:n_con) {
        results <- zeptosensPkg::get_target_score(
          wk = wk,
          wks = wks,
          dist_ind = dist_ind,
          inter = inter,
          n_dose = 1,
          n_prot = n_pro,
          proteomic_responses = proteomic_responses[i, ],
          max_dist = maxi_dist,
          n_perm = n_perm,
          cellLine = file_name,
          verbose = FALSE,
          fsFile = fs_value
        )
        ts[i, ] <- results$ts
        ts_p[i, ] <- results$pts
        ts_q[i, ] <- results$q
      }
      colnames(ts) <- colnames(drug_data)
      rownames(ts) <- rownames(drug_data)

      colnames(ts_p) <- colnames(drug_data)
      rownames(ts_p) <- rownames(drug_data)

      colnames(ts_q) <- colnames(drug_data)
      rownames(ts_q) <- rownames(drug_data)
    }

    if (ts_type == "pop") {

      # Calc Std
      drug_data[is.na(drug_data)] <- 0
      stdev <- zeptosensPkg::samp_sdev(nSample = nrow(drug_data), n_prot = ncol(drug_data), n_dose = 1, nX = drug_data)
      # normalization
      proteomic_responses <- drug_data
      for (i in 1:n_prot) {
        for (j in 1:nrow(proteomic_responses)) {
          proteomic_responses[j, i] <- (drug_data[j, i] / stdev[i])
        }
      }

      n_perm <- 25 # FIXME

      results <- zeptosensPkg::get_target_score(
        wk = wk,
        wks = wks,
        dist_ind = dist_ind,
        inter = inter,
        n_dose = 1,
        n_prot = n_pro,
        proteomic_responses = proteomic_responses,
        max_dist = maxi_dist,
        n_perm = n_perm,
        cellLine = file_name,
        verbose = FALSE,
        fsFile = fs_value
      )
      ts <- results$ts
      ts_p <- results$pts
      ts_q <- results$q

      colnames(ts) <- colnames(drug_data)
      colnames(ts_p) <- colnames(drug_data)
      colnames(ts_q) <- colnames(drug_data)
    }
    ts_r <- list(ts = ts, ts_p = ts_p, ts_q = ts_q)
    return(ts_r)
  })

  # data heatmap
  output$heatmap <- renderPlot({
    drug_file <- drug_dat()

    max_dat <- max(as.matrix(drug_file))
    min_dat <- min(as.matrix(drug_file))
    bk <- c(seq(min_dat, -0.01, by = 0.01), seq(0, max_dat, by = 0.01))
    data <- as.matrix(drug_file)
    pheatmap(data,
      scale = "none",
      color = c(
        colorRampPalette(colors = c("navy", "white"))(length(seq(min_dat, -0.01, by = 0.01))),
        colorRampPalette(colors = c("white", "firebrick3"))(length(seq(0, max_dat, by = 0.01)))
      ),
      legend_breaks = seq(min_dat, max_dat, 2), cellwidth = 2, cellheight = 2, fontsize = 2, fontsize_row = 2,
      breaks = bk
    )
  })

  #### DATA TABLE MODULE ----
  output$fs_value <- renderDataTable({
    data <- fs_dat()
    return(data)
  })

  output$anti_map <- renderDataTable({
    data <- anti_dat()
    return(data)
  })

  output$edgelist <- renderDataTable({
    data <- network_inferred()
    edgelist <- zeptosensPkg::create_sif_from_matrix(t.net = data$wk, genelist = colnames(data$wk))
    return(edgelist)
  })

  #### TEST MODULE ----
  output$test <- renderDataTable({
    data <- fs_value()
    return(data)
  })

  #### NETWORK VISUALIZATION MODULE ----
  #### TS HEATMAP MODULE ----
  output$tsheat <- renderPlot({
    ts_r <- ts_r()
    ts <- ts_r$ts
    max_dat <- max(as.matrix(ts))
    min_dat <- min(as.matrix(ts))
    bk <- c(seq(min_dat, -0.01, by = 0.01), seq(0, max_dat, by = 0.01))
    data <- as.matrix(ts)
    pheatmap(data,
      scale = "none",
      color = c(
        colorRampPalette(colors = c("navy", "white"))(length(seq(min_dat, -0.01, by = 0.01))),
        colorRampPalette(colors = c("white", "firebrick3"))(length(seq(0, max_dat, by = 0.01)))
      ),
      legend_breaks = seq(min_dat, max_dat, 2), cellwidth = 2, cellheight = 2, fontsize = 2, fontsize_row = 2,
      breaks = bk
    )
  })

  #### TS VOLCANO PLOT MODULE ----
  output$volcanoplot <- renderPlot({
    ts_r <- ts_r()
    n_line <- nline()
    ts <- ts_r$ts[n_line, ]
    ts_q <- ts_r$ts_q[nline, ]
    ts <- as.matrix(ts)
    p_adj <- as.matrix(ts_q)
    #
    if (nrow(p_adj) != nrow(ts)) {
      stop("ERROR: Tag of ts and q_value does not match.")
    }
    get_volcano_plot(ts = ts, q_value = ts_q, filename = rownames(ts_r)[n_line], path = "")
  })
}

shinyApp(server = server, ui = ui)
