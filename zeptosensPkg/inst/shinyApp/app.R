# Update Module
# library(rsconnect)
# rsconnect::setAccountInfo(name='hithut', token='B512E4272720432CD9C6F5E47388416D',
#   secret='/Y5XUfYMz02GGocWAcE/vQ+AsSrhpfEA8mBSk8Gl')
# rsconnect::deployApp()

library(shiny)
library(zeptosensPkg)

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

        actionButton("submit", label = "Submit", icon = NULL, width = NULL)
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
            "TargetScore",
            # heatmap of the data
            plotOutput("tsheat", height = 200, width = 1000)
          ),
          tabPanel(
            "TargetScore Volcano Plot",
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

  # Read in input files

  # Drug File
  results <- eventReactive(input$submit, {
    # Drug Data
    drug_file <- input$drug_data
    if (is.null(drug_file)) {
      return(NULL)
    }
    drug_dat <- read.csv(drug_file$datapath, header = TRUE)

    # Size of drug dat
    n_prot <- ncol(drug_dat)
    n_cond <- nrow(drug_dat)

    # Antibody Map
    antibody_map_file <- input$antibody
    if (is.null(antibody_map_file)) {
      anti_dat <- system.file("targetScoreData", "antibodyMapfile_08092019.txt", package = "zeptosensPkg")
    }
    if (!is.null(antibody_map_file)) {
      anti_dat <- read.csv(antibody_map_file$datapath, header = TRUE, stringsAsFactors = FALSE)
    }

    # FS File
    fs_file <- input$fs_file
    if (is.null(fs_file)) {
      return(NULL)
    }
    fs_dat <- read.csv(fs_file$datapath, header = TRUE, stringsAsFactors = FALSE)

    if (is.null(fs_dat)) {
      fs_value <- zeptosensPkg::get_fs_vals(
        n_prot = n_prot,
        proteomic_responses = drug_data
      )
    }
    if (!is.null(fs_dat)) {
      fs_value <- zeptosensPkg::get_fs_vals(
        n_prot = n_prot,
        proteomic_responses = drug_data,
        fs_value_file = fs_dat
      )
    }

    # ts_calc_type
    ts_type <- input$ts_calc_type

    # File name
    file_name <- input$filename

    # Line number
    n_line <- input$line

    # Max distance
    max_dist <- input$max_dist

    # Network algorithm
    network_algorithm <- input$network_algorithm

    # Proteomics dataset for network inference
    sig_file <- input$sig
    if (is.null(sig_file)) {
      return(NULL)
    }
    sig_dat <- read.csv(sig_file$datapath, header = TRUE, stringsAsFactors = FALSE)

    # choosing the way to construct reference network
    if (network_algorithm == "Bio") {
      # reference network
      network <- zeptosensPkg::predict_bio_network(
        n_prot = n_prot,
        proteomic_responses = drug_data,
        mab_to_genes = anti_data,
        max_dist = max_dist
      )
      wk <- network$wk
      wks <- network$wks
      dist_ind <- network$dist_ind
      inter <- network$inter
    }

    if (network_algorithm == "Dat") {
      network <- zeptosensPkg::predict_dat_network(
        data = sig_dat,
        n_prot = n_prot,
        proteomic_responses = drug_data,
        max_dist = max_dist
      )
      wk <- network$wk
      wks <- network$wks
      dist_ind <- network$dist_ind
      inter <- network$inter
    }

    if (network_algorithm == "Hyb") {
      # prior
      wk <- zeptosensPkg::predict_bio_network(
        n_prot = n_prot,
        proteomic_responses = drug_data,
        max_dist = max_dist,
        mab_to_genes = anti_data
      )$wk
      # Hyb
      network <- zeptosensPkg::predict_hybrid_network(
        data = sig_dat,
        prior = wk,
        n_prot = n_prot,
        proteomic_responses = drug_data
      )

      wk <- network$wk
      wks <- network$wks
      dist_ind <- network$dist_ind
      inter <- network$inter
    }
    network_inferred <- list(wk = wk, wks = wks, dist_ind = dist_ind, inter = inter)

    # Network inferred
    wk <- network_inferred$wk
    wks <- network_inferred$wks
    dist_ind <- network_inferred$dist_ind
    inter <- network_inferred$inter

    # Calculate Targetscore
    if (ts_type == "lnl") {
      drug_data[is.na(drug_data)] <- 0 # FIXME

      # Calc Std (Normalization request)
      # FIXME REMOVE?
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
      ts <- array(0, dim = c(n_cond, n_prot))
      ts_p <- array(0, dim = c(n_cond, n_prot))
      ts_q <- array(0, dim = c(n_cond, n_prot))

      for (i in 1:n_cond) {
        results <- zeptosensPkg::get_target_score(
          wk = wk,
          wks = wks,
          dist_ind = dist_ind,
          inter = inter,
          n_dose = 1,
          n_prot = n_prot,
          proteomic_responses = proteomic_responses[i, ],
          max_dist = max_dist,
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
        n_prot = n_prot,
        proteomic_responses = proteomic_responses,
        max_dist = max_dist,
        n_perm = n_perm,
        cellLine = file_name,
        verbose = FALSE,
        fs_file = fs_value
      )
      ts <- results$ts
      ts_p <- results$pts
      ts_q <- results$q

      colnames(ts) <- colnames(drug_data)
      colnames(ts_p) <- colnames(drug_data)
      colnames(ts_q) <- colnames(drug_data)
    }
    ts_r <- list(ts = ts, ts_p = ts_p, ts_q = ts_q)

    results <- list(
      ts_r = ts_r,
      drug_dat = drug_dat,
      fs_dat = fs_dat,
      anti_dat = anti_dat,
      fs_value = fs_value,
      network_inferred = network_inferred,
      n_line = n_line
    )
  })

  # OUTPUT ----

  ## Output heatmap
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

  # DATA TABLE MODULE ----
  output$fs_value <- renderDataTable({
    results <- results()
    return(results$fs_value)
  })

  output$anti_map <- renderDataTable({
    results <- results()
    return(results$anti_dat)
  })

  output$edgelist <- renderDataTable({
    results <- results()
    network_inferred <- results$network_inferred
    edgelist <- zeptosensPkg::create_sif_from_matrix(
      t.net = network_inferred$wk,
      genelist = colnames(network_inferred$wk)
    )
    return(edgelist)
  })

  #### TEST MODULE ----
  output$test <- renderDataTable({
    results <- results()
    return(results$fs_value)
  })

  #### NETWORK VISUALIZATION MODULE ----
  #### TS HEATMAP MODULE ----
  output$tsheat <- renderPlot({
    results <- results()
    ts_r <- results$ts_r
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
    results <- results()
    ts_r <- results$ts_r
    n_line <- results$n_line
    ts <- ts_r$ts[n_line, ]
    ts_q <- ts_r$ts_q[n_line, ]
    ts <- as.matrix(ts)
    p_adj <- as.matrix(ts_q)

    if (nrow(p_adj) != nrow(ts)) {
      stop("ERROR: Number of rows in TS and adjusted p-value do not match.")
    }

    get_volcano_plot(ts = ts, q_value = ts_q, filename = rownames(ts_r)[n_line], path = "")
  })
}

shinyApp(server = server, ui = ui)
