# install neccessary libraries

#if(!require("devtools")) { install.packages("devtools") }
#if(!require("morpheus")) { devtools::install_github('cmap/morpheus.R') }
#if(!require("zeptosensPkg")) { devtools::install_github("HepingWang/zeptosensPkg/zeptosensPkg") }

library(shiny)
library(markdown)
library(DT)
library(zeptosensPkg)
library(pheatmap)
library(morpheus)
library(plotly)
library(devtools)
source("modal.R")


# UI ----
ui <- navbarPage(
  "Target Score",
  header = list(
    # Add/run startup Javascript
    tags$head(tags$script(onload_js)),
    # Use JQuery (built into Shiny) calls to show/hide modal based on message
    tags$head(includeScript("www/js/showLoading.js"))
  ),
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
        fileInput("ts_result_file", "TargetScore Result File (.rds; Visualization Only)",
                  buttonLabel = "Browse...",
                  placeholder = "No file selected",
                  accept = ".rds"
        ),
        fileInput("drug_data_file", "Perturbation Response File (.csv; REQUIRED)",
          buttonLabel = "Browse...",
          placeholder = "No file selected",
          accept = ".csv"
        ),
        fileInput("antibody", "Mapping File (.csv or blank)",
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
          label = "Target Score Calculation Type:",
          choices = c(
            "Line by Line" = "line_by_line",
            "Pooled" = "pooled"
          ),
          selected = "line_by_line"
        ),

        # Number of permutations
        numericInput("n_perm", "Permutation Number", value = 25, min = 25, max = 2000), # FIXME

        # max distance of protein network
        numericInput("max_dist", "Maximum Network Distance", "1"),

        actionButton("submit", label = "Submit", icon = NULL, width = NULL),

        hr(),
        helpText("Example: ", a(href="data/sample_breast_cancer_perturbation_data.zip", target="_blank", download="sample_breast_cancer_perturbation_data.zip", "Sample Breast Dataset (.zip)"))
        # helpText(
        #   a("Perturbation Response Example", href = "data/BT474.csv", target = "_blank"),
        #   br(),
        #   a("Mapping Example", href = "data/antibodyMapfile.txt", target = "_blank"),
        #   br(),
        #   a("Background Network Example", href = "data/TCGA_BRCA_L4.csv", target = "_blank"),
        #   br(),
        #   a("Functional Score Example", href = "data/Cosmic.txt", target = "_blank"),
        #   br(),
        #   a("Pre-computed Results for Visualization Example", href = "data/BT474_results_example.rds", target = "_blank"),
        # ),
      ),
      # Results showing
      mainPanel(
        # Results showing in Tabs (can use navlistPanel to show on left)
        tabsetPanel(
          tabPanel(
            "Protein Mapping/Functional Scores",
             tableOutput("fs_dat")
            #DT::dataTableOutput("fs_dat")
          ),
          tabPanel(
            "Input Data Heatmap",
            # heatmap of the data
            morpheusOutput("heatmap")
          ),
          tabPanel(
            "Network Edgelist",
            # edgelist of the data
            tableOutput("edgelist")
            #DT::dataTableOutput("edgelist")
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
            numericInput("condition_number", "Condition (Input Row Number)", value = 1, min = 1, step = 1), # FIXME
            br(),
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
      includeMarkdown("www/ts_about.md"),
      h1("Version"),
      p(paste0("TargetScore: ", packageVersion("zeptosensPkg")))
    )
  ),
  loading_modal("Calculating TargetScore ...")
)

# SERVER ----
server <- function(input, output, session) {
  results <- eventReactive(input$submit, {
    cat("DEBUG\n")
    
    #Result Load in or Calculate start
    if(!is.null(input$ts_result_file)){
      ts_result_file <- input$ts_result_file
      results <- readRDS(ts_result_file$datapath)
    } else if(!is.null(input$drug_data_file)) {
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
        # fs_file <- system.file("targetScoreData", "fs.csv", package = "zeptosensPkg")
        fs_data <- read.csv(fs_file$datapath, header = TRUE, stringsAsFactors = FALSE)
        
        fs_dat <- zeptosensPkg::get_fs_vals(
          n_prot = n_prot,
          proteomic_responses = proteomic_responses,
          mab_to_genes = mab_to_genes,
          fs_override = fs_data
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
      
      # DEBUG
      cat("D1X: ", as.character(Sys.time()), "\n")
      
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
        validate(
          need(is.null(sig_file), "ERROR: REQUIRED: Background Network File for Hybrid and 
           Graphical LASSO Network Construction Algorithms")
        )
        
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
      
      # DEBUG
      cat("D1Y: ", as.character(Sys.time()), "\n")
      
      # Calculate Targetscore
      if (ts_type == "line_by_line") {
        proteomic_responses[is.na(proteomic_responses)] <- 0 # FIXME
        
        # Calc Std (Normalization request)
        # FIXME REMOVE? there should be no na's included in preoteomic responses. What if so?
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
          n_dose = nrow(proteomic_responses),
          n_prot = n_prot,
          proteomic_responses = proteomic_responses,
          n_perm = n_perm,
          verbose = FALSE,
          fs_dat = fs_dat
        )
        ts <- results$ts
        ts_p <- results$pts
        ts_q <- results$q
        
        names(ts) <- colnames(proteomic_responses)
        names(ts_p) <- colnames(proteomic_responses)
        names(ts_q) <- colnames(proteomic_responses)
      }
      ts_r <- list(ts = ts, ts_p = ts_p, ts_q = ts_q)
      
      # DEBUG
      cat("D2: ", str(ts_r), "\n")
      
      results <- list(
        ts_r = ts_r,
        proteomic_responses = proteomic_responses,
        fs_dat = fs_dat,
        mab_to_genes = mab_to_genes,
        network = network
      )
      
      # DEBUG
      saveRDS(results, "shiny_results_example.rds")
    } else {
      validate(
        need((!is.null(input$drug_data_file) | !is.null(input$ts_result_file)), "ERROR: REQUIRED: Drug Response File or TargetScore Result File")
      )
    }
    
    return(results)
  })

  # OBSERVERS ----
  # Start showing loading marker
  observeEvent(input$submit, {
    cat("SUBMITTED\n")
    session$sendCustomMessage(type = "showLoading", list(show = TRUE))
  })

  # Observe reactive variable and send message to Javascript code
  observe({
    results <- results()
    if (length(results) > 0) {
      cat("FINISHED\n")
      session$sendCustomMessage(type = "showLoading", list(show = FALSE))
    }
  })

  # OUTPUT ----

  ## Output heatmap
  output$heatmap <- renderMorpheus({
    results <- results()
    proteomic_responses <- results$proteomic_responses

    drug_dat_mat <- as.matrix(proteomic_responses)

    morpheus::morpheus(drug_dat_mat, dendrogram = "none")
  })

  # DATA TABLE MODULE ----
   output$fs_dat <- renderTable({
  #output$fs_dat <- DT::renderDataTable({
    results <- results()
    # return(results$fs_dat)

    t1 <- results$fs_dat
    t2 <- results$mab_to_genes
    dat <- merge(t1, t2, by.x = "prot", by.y = "AntibodyLabel")
    
    return(dat)
  })

  output$edgelist <- renderTable({
  #output$edgelist <- DT::renderDataTable({
    results <- results()
    network <- results$network
    edgelist <- zeptosensPkg::create_sif_from_matrix(
      t_net = network$wk,
      col_genelist = colnames(network$wk),
      row_genelist = rownames(network$wk)
    )
    return(edgelist)
  })

  #### TEST MODULE ----
  output$test <- renderTable({
  #output$test <- DT::renderDataTable({
    results <- results()
    tmp <- data.frame(a = 1, b = 2, stringsAsFactors = FALSE)
    return(tmp)
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

    morpheus::morpheus(ts_mat, dendrogram = "none")
  })

  #### TS VOLCANO PLOT MODULE ----
  # output$volcano_plot <- renderPlot({
  output$volcano_plot <- renderPlotly({
    # Line number
    condition_number <- input$condition_number
    
    
    calc_type <- input$ts_calc_type
    results <- results()
    ts_r <- results$ts_r
    
    validate(need(condition_number <= nrow(results$ts_r$ts), paste0("Total Available Conditions in Output: ", nrow(results$ts_r$ts))))
    
    if (calc_type == "line_by_line") {
      ts <- ts_r$ts[condition_number, ]
      ts_q <- ts_r$ts_q[condition_number, ]
    } else if (calc_type == "pooled") {
      ts <- ts_r$ts
      ts_q <- ts_r$ts_q
    }
    ts <- as.matrix(ts)
    p_adj <- as.matrix(ts_q)

    if (nrow(p_adj) != nrow(ts)) {
      stop("ERROR: Number of rows in TS and adjusted p-value do not match.")
    }

    p1 <- get_volcano_plot(
      ts = ts, q_value = ts_q, filename = rownames(ts_r)[condition_number], path = "",
      include_labels = FALSE, save_output = FALSE
    )
    
    # Add title
    p1 <- p1 + ggtitle(rownames(results$ts_r$ts)[condition_number])
    
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
