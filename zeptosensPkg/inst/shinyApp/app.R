# Update Module
# library(rsconnect)
# rsconnect::setAccountInfo(name='hithut', token='B512E4272720432CD9C6F5E47388416D',
#   secret='/Y5XUfYMz02GGocWAcE/vQ+AsSrhpfEA8mBSk8Gl')
# rsconnect::deployApp()

library(shiny)
library(zeptosensPkg)
library(pheatmap)
library(morpheus)
library(plotly)

# UI ----
ui <- navbarPage(
  "Target Score",
  tabPanel(
    "Overview",
    mainPanel(
      includeMarkdown("www/ts_intro_p1.md")
      # FIXME ADD INTRO.JPG: 1) redraw horizontal (preferred)
      # or 2) https://stackoverflow.com/questions/31603577/two-column-layout-with-markdown
    )
  ),
  tabPanel(
    "Run",
    sidebarLayout(
      sidebarPanel(
        width = 4,
        fileInput("antibody", "Antibody Mapping File (.csv or blank)",
          buttonLabel = "Browse...",
          placeholder = "No file selected",
          accept = ".csv"
        ),
        fileInput("sig", "Background Network File (.csv or blank)",
          buttonLabel = "Browse...",
          placeholder = "No file selected",
          accept = ".csv"
        ),
        fileInput("fs_file", "Functional Score File (.csv or blank)",
          buttonLabel = "Browse...",
          placeholder = "No file selected",
          accept = ".csv"
        ),
        fileInput("drug_data_file", "Drug Perturbation Response File (.csv)",
          buttonLabel = "Browse...",
          placeholder = "No file selected",
          accept = ".csv"
        ),
        selectInput("network_algorithm",
          label = "Network Construction Algorithms:",
          choices = c(
            "Literature-Based" = "bio",
            "Hybrid" = "hybrid",
            "Graphical LASSO" = "dat"
          ),
          selected = "bio"
        ),
        selectInput("ts_calc_type",
          label = "Target Score Calculation type:",
          choices = c(
            "Line by Line" = "line_by_line",
            "Pooled" = "pooled"
          ),
          selected = "line_by_line"
        ),

        # Number of permutations
        numericInput("n_perm", "Number of permutations", value = 25, min = 25, max = 2000), # FIXME

        # max distance of protein network
        numericInput("max_dist", "Maximum Protein Distance", "1"),

        actionButton("submit", label = "Submit", icon = NULL, width = NULL)
      ),
      # Results showing
      mainPanel(
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
            dataTableOutput("fs_dat")
          ),
          tabPanel(
            "Input Data Heatmap",
            # heatmap of the data
            plotOutput("heatmap", height = 200, width = 1000)
          ),
          tabPanel(
            "Network Edgelist",
            # edgelist of the data
            dataTableOutput("edgelist")
          ),
          tabPanel(
            "TargetScore Heatmap",
            # heatmap of the data
            # plotOutput("tsheat", height = 200, width = 1000)
            morpheusOutput("tsheat_morpheus")
          ),
          tabPanel(
            "TargetScore Volcano Plot",
            # volcano plot line choice
            numericInput("line_number", "Line Number", value = 1, min = 1), # FIXME
            plotlyOutput("volcano_plot")
          )
        )
      )
    )
  ),
  tabPanel(
    "Help",
    mainPanel(
      includeMarkdown("www/ts_intro_p2.md"),
      includeMarkdown("www/ts_input.md")
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
  results <- eventReactive(input$submit, {
    cat("DEBUG\n")

    validate(
      need(input$drug_data_file, "ERROR: REQUIRED: Drug Response File")
    )

    # Drug Data
    drug_data_file <- input$drug_data_file
    drug_data_file <- drug_data_file$datapath
    # DEBUG
    # drug_data_file <- system.file("test_data", "BT474.csv", package = "zeptosensPkg")
    cat("X: ", drug_data_file, "\n")

    # Read drug dataset, NOTE: must have row names
    proteomic_responses <- read.csv(drug_data_file, row.names = 1, header = TRUE)

    # DEBUG
    cat("D1: ", str(proteomic_responses), "\n")

    # Size of drug dat
    n_prot <- ncol(proteomic_responses)
    n_cond <- nrow(proteomic_responses)

    # Antibody Map
    antibody_map_file <- input$antibody
    if (is.null(antibody_map_file)) {
      antibody_map_file <- system.file("targetScoreData", "antibodyMapfile_08092019.txt", package = "zeptosensPkg")
      mab_to_genes <- read.table(antibody_map_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    } else {
      antibody_map_file <- antibody_map_file$datapath
      mab_to_genes <- read.csv(antibody_map_file, header = TRUE, stringsAsFactors = FALSE)
    }

    # FS File
    fs_file <- input$fs_file
    if (!is.null(fs_file)) {
      fs_file <- system.file("targetScoreData", "fs.csv", package = "zeptosensPkg")
      fs_dat <- read.csv(fs_file, header = TRUE, stringsAsFactors = FALSE)

      fs_dat <- zeptosensPkg::get_fs_vals(
        n_prot = n_prot,
        proteomic_responses = proteomic_responses,
        mab_to_genes = mab_to_genes,
        fs_dat = fs_dat
      )
    } else {
      fs_dat <- zeptosensPkg::get_fs_vals(
        n_prot = n_prot,
        proteomic_responses = proteomic_responses,
        mab_to_genes = mab_to_genes
      )
    }

    # ts_calc_type
    ts_type <- input$ts_calc_type

    # Max distance
    max_dist <- input$max_dist

    # Network algorithm
    network_algorithm <- input$network_algorithm

    # Number of permutations
    n_perm <- input$n_perm

    # choosing the way to construct reference network
    if (network_algorithm == "bio") {
      # reference network
      network <- zeptosensPkg::predict_bio_network(
        n_prot = n_prot,
        proteomic_responses = proteomic_responses,
        mab_to_genes = mab_to_genes,
        max_dist = max_dist
      )
      wk <- network$wk
      wks <- network$wks
      dist_ind <- network$dist_ind
      inter <- network$inter
    }

    if (network_algorithm == "dat" || network_algorithm == "hybrid") {
      # Proteomics dataset for network inference
      sig_file <- input$sig
      if (is.null(sig_file)) {
        return(NULL)
      }
      # FIXME stringsAsFactors for anything else?
      sig_dat <- read.csv(sig_file$datapath, header = TRUE, stringsAsFactors = FALSE)

      if (network_algorithm == "dat") {
        network <- zeptosensPkg::predict_dat_network(
          data = sig_dat,
          n_prot = n_prot,
          proteomic_responses = proteomic_responses,
          max_dist = max_dist
        )
        wk <- network$wk
        wks <- network$wks
        dist_ind <- network$dist_ind
        inter <- network$inter
      }

      if (network_algorithm == "hybrid") {
        # prior
        wk <- zeptosensPkg::predict_bio_network(
          n_prot = n_prot,
          proteomic_responses = proteomic_responses,
          max_dist = max_dist,
          mab_to_genes = mab_to_genes
        )$wk
        # Hyb
        network <- zeptosensPkg::predict_hybrid_network(
          data = sig_dat,
          prior = wk,
          n_prot = n_prot,
          proteomic_responses = proteomic_responses
        )

        wk <- network$wk
        wks <- network$wks
        dist_ind <- network$dist_ind
        inter <- network$inter
      }
    }

    # Calculate Targetscore
    if (ts_type == "line_by_line") {
      proteomic_responses[is.na(proteomic_responses)] <- 0 # FIXME

      # Calc Std (Normalization request)
      # FIXME REMOVE?
      # stdev <- zeptosensPkg::samp_sdev(nSample=nrow(proteomic_responses),
      #   n_prot=ncol(proteomic_responses),n_dose=1,nX=proteomic_responses)
      # #normalization
      # proteomic_responses<- proteomic_responses
      # for(i in 1:n_prot) {
      #   for (j in 1:nrow(proteomic_responses)) {
      #     proteomic_responses[j,i] <- (proteomic_responses[j,i]/stdev[i])
      #   }
      # }
      #
      proteomic_responses <- proteomic_responses
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
          n_perm = n_perm,
          verbose = FALSE,
          fs_dat = fs_dat
        )
        ts[i, ] <- results$ts
        ts_p[i, ] <- results$pts
        ts_q[i, ] <- results$q
      }
      colnames(ts) <- colnames(proteomic_responses)
      rownames(ts) <- rownames(proteomic_responses)

      colnames(ts_p) <- colnames(proteomic_responses)
      rownames(ts_p) <- rownames(proteomic_responses)

      colnames(ts_q) <- colnames(proteomic_responses)
      rownames(ts_q) <- rownames(proteomic_responses)
    }

    if (ts_type == "pooled") {
      # Calc Std
      stdev <- zeptosensPkg::samp_sdev(
        n_sample = nrow(proteomic_responses),
        n_prot = ncol(proteomic_responses),
        n_dose = 1,
        n_x = proteomic_responses,
        replace_missing = TRUE
      )
      # normalization
      proteomic_responses <- proteomic_responses
      for (i in 1:n_prot) {
        for (j in 1:nrow(proteomic_responses)) {
          proteomic_responses[j, i] <- (proteomic_responses[j, i] / stdev[i])
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
        n_perm = n_perm,
        verbose = FALSE,
        fs_dat = fs_dat
      )
      ts <- results$ts
      ts_p <- results$pts
      ts_q <- results$q

      colnames(ts) <- colnames(proteomic_responses)
      colnames(ts_p) <- colnames(proteomic_responses)
      colnames(ts_q) <- colnames(proteomic_responses)
    }
    ts_r <- list(ts = ts, ts_p = ts_p, ts_q = ts_q)

    # DEBUG
    cat("D2: ", str(ts_r), "\n")

    results <- list(
      ts_r = ts_r,
      proteomic_responses = proteomic_responses,
      fs_dat = fs_dat,
      mab_to_genes = mab_to_genes,
      fs_dat = fs_dat,
      network = network
    )
  })

  # OUTPUT ----

  ## Output heatmap
  output$heatmap <- renderPlot({
    results <- results()
    proteomic_responses <- results$proteomic_responses

    drug_dat_mat <- as.matrix(proteomic_responses)

    max_dat <- max(drug_dat_mat)
    min_dat <- min(drug_dat_mat)
    bk <- c(seq(min_dat, -0.01, by = 0.01), seq(0, max_dat, by = 0.01))

    pheatmap(drug_dat_mat,
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
  output$fs_dat <- renderDataTable({
    results <- results()
    return(results$fs_dat)
  })

  output$anti_map <- renderDataTable({
    results <- results()
    return(results$mab_to_genes)
  })

  output$edgelist <- renderDataTable({
    results <- results()
    network <- results$network
    edgelist <- zeptosensPkg::create_sif_from_matrix(
      t_net = network$wk,
      genelist = colnames(network$wk)
    )
    return(edgelist)
  })

  #### TEST MODULE ----
  output$test <- renderDataTable({
    results <- results()
    return("CALCULATION DONE")
  })

  #### NETWORK VISUALIZATION MODULE ----
  #### TS HEATMAP MODULE ----
  output$tsheat <- renderPlot({
    results <- results()
    ts_r <- results$ts_r
    ts <- ts_r$ts
    ts_mat <- as.matrix(ts)

    max_dat <- max(ts_mat)
    min_dat <- min(ts_mat)
    bk <- c(seq(min_dat, -0.01, by = 0.01), seq(0, max_dat, by = 0.01))

    pheatmap(ts_mat,
      scale = "none",
      color = c(
        colorRampPalette(colors = c("navy", "white"))(length(seq(min_dat, -0.01, by = 0.01))),
        colorRampPalette(colors = c("white", "firebrick3"))(length(seq(0, max_dat, by = 0.01)))
      ),
      legend_breaks = seq(min_dat, max_dat, 2), cellwidth = 2, cellheight = 2, fontsize = 2, fontsize_row = 2,
      breaks = bk
    )
  })

  output$tsheat_morpheus <- renderMorpheus({
    results <- results()
    ts_r <- results$ts_r
    ts <- ts_r$ts
    ts_mat <- as.matrix(ts)

    max_dat <- max(ts_mat)
    min_dat <- min(ts_mat)
    bk <- c(seq(min_dat, -0.01, by = 0.01), seq(0, max_dat, by = 0.01))

    x <- matrix(rnorm(200), 20)
    y <- data.frame(a = letters[1:10], b = rep(c("g", "h"), 5), stringsAsFactors = FALSE)
    morpheus::morpheus(x, columnAnnotations = y)
  })

  #### TS VOLCANO PLOT MODULE ----
  # output$volcano_plot <- renderPlot({
  output$volcano_plot <- renderPlotly({
    # Line number
    line_number <- input$line_number

    results <- results()
    ts_r <- results$ts_r
    ts <- ts_r$ts[line_number, ]
    ts_q <- ts_r$ts_q[line_number, ]
    ts <- as.matrix(ts)
    p_adj <- as.matrix(ts_q)

    if (nrow(p_adj) != nrow(ts)) {
      stop("ERROR: Number of rows in TS and adjusted p-value do not match.")
    }

    p1 <- get_volcano_plot(
      ts = ts, q_value = ts_q, filename = rownames(ts_r)[line_number], path = "",
      include_labels = FALSE, save_output = FALSE
    )
    # g1 <- ggplotly(p1, width=plotWidth, height=plotHeight, tooltip=tooltipCol) # need tooltip
    g1 <- plotly::ggplotly(p1)
    # g1 <- layout(g1, margin=list(t = 75))
    g2 <- plotly::config(
      p = g1, cloud = FALSE, displaylogo = FALSE,
      modeBarButtonsToRemove = c(
        "select2d", "sendDataToCloud", "pan2d", "resetScale2d",
        "hoverClosestCartesian", "hoverCompareCartesian",
        "lasso2d", "zoomIn2d", "zoomOut2d"
      )
    )
    # p1
    g2
  })
}

shinyApp(server = server, ui = ui)
