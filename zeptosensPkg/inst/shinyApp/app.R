# Update Module
# library(rsconnect)
# rsconnect::setAccountInfo(name='hithut', token='B512E4272720432CD9C6F5E47388416D',
#   secret='/Y5XUfYMz02GGocWAcE/vQ+AsSrhpfEA8mBSk8Gl')
# rsconnect::deployApp()

library(pheatmap)
library(shiny)
library(glasso)
library(zeptosensPkg)
library(ggplot2)
library(ggrepel)
library(markdown)

n_perm <- 25

# UI ----
ui <- navbarPage(
  "Target Score",
  tabPanel(
    "Get Started",
    mainPanel(
      includeMarkdown("www/ts_intro_p1.md"),
      includeMarkdown("www/ts_intro_p2.md")
      # FIXME ADD INTRO.JPG: 1) redraw horizontal (preferred)
      # or 2) https://stackoverflow.com/questions/31603577/two-column-layout-with-markdown
    )
  ),
  tabPanel(
    "Input File Descriptions",
    mainPanel(
      includeMarkdown("www/ts_input.md")
    )
  ),
  tabPanel(
    "App",
    sidebarLayout(
      sidebarPanel(
        width = 3,
        fileInput("antibody", "Antibody Mapping File (.csv)",
          buttonLabel = "Browse...",
          placeholder = "No file selected",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        ),
        fileInput("sig", "Background Network File (.csv)",
          buttonLabel = "Browse...",
          placeholder = "No file selected",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        ),
        fileInput("fs_file", "Functional Score File (.csv)",
          buttonLabel = "Browse...",
          placeholder = "No file selected",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        ),
        fileInput("drug_data", "Drug Perturbation Response File (.csv)",
          buttonLabel = "Browse...",
          placeholder = "No file selected",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        ),
        selectInput("network_algorithm",
          label = "Network Construction Algorithms:",
          choices = c(
            "Hybrid-driven" = "Hyb",
            "Biology-inferral" = "Bio",
            "data-driven" = "Dat"
          ),
          selected = NULL
        ),
        selectInput("ts_calc_type",
          label = "Target Score Calculation type:",
          choices = c(
            "Line by Line" = "lnl",
            "Pooled" = "pop"
          ),
          selected = NULL
        ),
        # volcano plot line choice
        numericInput("line", "Line Number", "1"),
        # max distance of protein network
        numericInput("max_dist", "Maximum Protein Distance", "1"),

        actionButton("Submit", label = "Submit", icon = NULL, width = NULL)
      ),
      # Results showing
      mainPanel(
        tags$h2("Guideline:"),
        hr(),
        tags$a(href = "http://www.git.com", "TargetScore Package"),

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
            "Network Edgelist",
            # edgelist of the data
            dataTableOutput("edgelist")
          ),
          tabPanel(
            "Targetscore",
            # heatmap of the data
            plotOutput("tsheat", height = 200, width = 1000)
          ),
          tabPanel(
            "Targetscore Volcano Plot",
            # Plots
            plotOutput("volcanoplot")
          )
        )
      )
    )
  ),
  tabPanel(
    "About",
    mainPanel(
      includeMarkdown("www/ts_about.md")
    )
  )
)

# SERVER ----
server <- function(input, output, session) {

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
  # CalcType <- eventReactive(input$Submit,{input$ts_calc_type})

  # load in system file and Default value and the Default File

  # file from input

  # Drug File
  drug_dat <- reactive({
    drug_file <- input$drug_data
    if (is.null(drug_file)) {
      return(NULL)
    }
    drug_dat <- read.csv(drug_file$datapath, header = TRUE)
    return(drug_dat)
  })

  # AntibodyMap File (Default at System file)
  anti_dat <- reactive({
    antibody_map_file <- input$antibody
    if (is.null(antibody_map_file)) {
      anti_dat <- system.file("targetScoreData", "antibodyMapfile_08092019.txt", package = "zeptosensPkg")
    }
    if (!is.null(antibody_map_file)) {
      anti_dat <- read.csv(antibody_map_file$datapath, header = TRUE, stringsAsFactors = FALSE)
    }
    return(anti_dat)
  })

  # FS File
  fs_dat <- reactive({
    fs_file <- input$fs_file
    if (is.null(fs_file)) {
      return(NULL)
    }
    fs_dat <- read.csv(fs_file$datapath, header = TRUE, stringsAsFactors = FALSE)

    return(fs_dat)
  })

  # ts_calc_type
  ts_type <- reactive({
    ts_type <- input$ts_calc_type
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
    sig_dat <- read.csv(sig_file$datapath, header = TRUE, stringsAsFactors = FALSE)
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
    nline <- input$line
    return(nline)
  })

  # max_dist
  max_dist <- reactive({
    max_dist <- input$max_dist
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
        mab_to_genes = anti_data,
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
        mab_to_genes = anti_data
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
    n_pro <- n_prot()
    fs_dat <- fs_dat()
    if (is.null(fs_dat)) {
      fs_value <- zeptosensPkg::get_fs_vals(
        n_prot = n_pro,
        proteomic_responses = drug_data
      )
    }
    if (!is.null(fs_dat)) {
      fs_value <- zeptosensPkg::get_fs_vals(
        n_prot = n_pro,
        proteomic_responses = drug_data,
        fs_value_file = fs_dat
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
    fs_value <- fs_dat()
    file_name <- file_name()
    maxi_dist <- max_dist()

    # Network inferred
    wk <- network_inferred$wk
    wks <- network_inferred$wks
    dist_ind <- network_inferred$dist_ind
    inter <- network_inferred$inter

    if (ts_type == "lnl") {
      drug_data[is.na(drug_data)] <- 0 # FIXME

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
          cell_line = file_name,
          verbose = FALSE,
          fs_file = fs_value
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
      stdev <- zeptosensPkg::samp_sdev(
        n_sample = nrow(drug_data),
        n_prot = ncol(drug_data),
        n_dose = 1,
        n_x = drug_data,
        replace_missing = TRUE
      )
      # normalization
      proteomic_responses <- drug_data
      for (i in 1:n_prot) {
        for (j in 1:nrow(proteomic_responses)) {
          proteomic_responses[j, i] <- (drug_data[j, i] / stdev[i])
        }
      }

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

    if (nrow(p_adj) != nrow(ts)) {
      stop("ERROR: Number of rows in TS and adjusted p-value do not match.")
    }

    get_volcano_plot(ts = ts, q_value = ts_q, filename = rownames(ts_r)[n_line], path = "")
  })
}

shinyApp(server = server, ui = ui)
