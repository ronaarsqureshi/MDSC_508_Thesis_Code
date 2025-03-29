# Load libraries
library(shiny)
library(readxl)
library(reactable)
library(openxlsx)   # For exporting Excel files
library(plotly)     # For interactive plots
library(ggtree)     # For phylogenetic trees
library(ape)        # For tree structures
library(DECIPHER)   # For alignment and distance calculation
library(Biostrings) # For reading fasta files and translation
library(dplyr)      # For data manipulation
library(msa)        # For multiple sequence alignment
library(phangorn)   # For MLE phylogenetic tree estimation
library(leaflet)    # For geospatial mapping
library(tidyr)      # For pivot_longer

# Load custom filtering function 
source("/Users/ronaarqureshi/Desktop/asv_filter.R")---------------------------------------------------------------------(ADJUST PATH HERE IF NEED)

ui <- fluidPage(
  titlePanel("ASV Metadata Filter"),
  
  # Button that appears only after extra tabs are inserted
  uiOutput("return_button_ui"),
  
  tabsetPanel(
    id = "main_tabs",
    
    # 1) Table View
    tabPanel(
      "Table View",
      sidebarLayout(
        sidebarPanel(
          fileInput("metadata_file", "Upload Metadata File", accept = c(".xlsx")),
          
          # Single species filter for the entire app
          selectInput("table_species_filter", "Filter by Species:", choices = NULL, selected = ""),
          
          numericInput("min_reads", "Minimum Read Number:", value = 10, min = 0),
          numericInput("min_percentage", "Minimum Read Percentage:", value = 0.1, min = 0, max = 100),
          numericInput("min_samples", "Minimum Sample Count:", value = 5, min = 0),
          downloadButton("download_filtered", "Download Filtered Data"),
          actionButton("apply_filters", "Apply Filters")
        ),
        mainPanel(
          reactableOutput("asv_table")
        )
      )
    ),
    
    # 2)Allele Frequency Plot
    tabPanel(
      "Allele Frequency Plot",
      mainPanel(
        div(
          style = "overflow-x: scroll;",
          plotlyOutput("allele_freq_plot", height = "600px", width = "1200px")
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  metadata <- reactive({
    if (is.null(input$metadata_file)) {
      list(
        sample_md = read_excel("/Users/ronaarqureshi/Desktop/metadata.xlsx", sheet = "SAMPLE_MD"),
        asv_md    = read_excel("/Users/ronaarqureshi/Desktop/metadata.xlsx", sheet = "ASV_MD"),
        asv_read  = read_excel("/Users/ronaarqureshi/Desktop/metadata.xlsx", sheet = "ASV_READ")
      )
    } else {
      list(
        sample_md = read_excel(input$metadata_file$datapath, sheet = "SAMPLE_MD"),
        asv_md    = read_excel(input$metadata_file$datapath, sheet = "ASV_MD"),
        asv_read  = read_excel(input$metadata_file$datapath, sheet = "ASV_READ")
      )
    }
  })
  
  #Populate species filter choices (for the entire app)
  observe({
    req(metadata()$asv_md)
    updateSelectInput(
      session, "table_species_filter",
      choices  = c("", unique(metadata()$asv_md$ASV_SPECIES)),
      selected = ""
    )
  })
  
  filtered_data <- reactive({
    fd <- asv_filter(
      metadata()$asv_read,
      metadata()$asv_md,
      metadata()$sample_md,
      min_read_num   = input$min_reads,
      min_read_pct   = input$min_percentage,
      min_sample_num = input$min_samples
    )
    
    # Apply species filter if selected
    if (!is.null(input$table_species_filter) && input$table_species_filter != "") {
      fd$md <- fd$md %>% filter(ASV_SPECIES == input$table_species_filter)
    }
    fd
  })
  
  # Show filtered table
  output$asv_table <- renderReactable({
    req(filtered_data()$md)
    explanations <- list(
      ASV_ID      = "Unique identifier for the ASV",
      ASV_SPECIES = "Species name",
      READ_NUM    = "Number of reads",
      SAMPLE_ID   = "Sample identifier"
    )
    col_defs <- lapply(names(filtered_data()$md), function(col) {
      if (col %in% names(explanations)) {
        colDef(header = function(value) tags$span(title = explanations[[col]], value))
      } else {
        colDef()
      }
    })
    names(col_defs) <- names(filtered_data()$md)
    
    reactable(filtered_data()$md, pagination = TRUE, searchable = TRUE, columns = col_defs)
  })
  
  # Download filtered data
  output$download_filtered <- downloadHandler(
    filename = function() { "metadata_filtered.xlsx" },
    content  = function(file) {
      write.xlsx(filtered_data(), file)
    }
  )
  
  
  output$allele_freq_plot <- renderPlotly({
    req(input$table_species_filter)
    
    # Use filtered_data()$md for the selected species
    asv_data <- filtered_data()$md
    if (nrow(asv_data) == 0) {
      return(plotly_empty(type = "scatter", mode = "markers"))
    }
    
    selected_species <- input$table_species_filter
    asv_data <- asv_data %>% filter(ASV_SPECIES == selected_species)
    
    if (nrow(asv_data) == 0) {
      return(plotly_empty(type = "scatter", mode = "markers"))
    }
    
    # Use the original (unfiltered) asv_read for pivot, but only keep ASVs from asv_data
    asv_read_species <- metadata()$asv_read %>%
      pivot_longer(cols = -SAMPLE_ID, names_to = "ASV_ID", values_to = "READ_NUM") %>%
      filter(ASV_ID %in% asv_data$ASV_ID, READ_NUM > 0)
    
    if (nrow(asv_read_species) == 0) {
      return(plotly_empty(type = "scatter", mode = "markers"))
    }
    
    freq_data <- asv_read_species %>%
      group_by(SAMPLE_ID) %>%
      mutate(Total_READS = sum(READ_NUM)) %>%
      ungroup() %>%
      mutate(Frequency = (READ_NUM / Total_READS) * 100) %>%
      select(SAMPLE_ID, ASV_ID, Frequency)
    
    plot_ly(freq_data, x = ~SAMPLE_ID, y = ~Frequency, color = ~ASV_ID, type = 'bar') %>%
      layout(
        title  = paste("Allele Frequency per Sample for", selected_species),
        xaxis  = list(title = "Sample ID"),
        yaxis  = list(title = "Frequency (%)"),
        barmode= 'stack'
      )
  })
  
  
  observe({
    req(filtered_data()$md)
    species_asvs <- filtered_data()$md %>%
      pull(ASV_ID) %>%
      unique()
    
    updateSelectInput(session, "map_asv_select",
                      choices = species_asvs,
                      selected = species_asvs)
  })
  
  
  tabsInserted <- reactiveVal(FALSE)
  
  # new_dataset is created after "Apply Filters"
  new_dataset <- eventReactive(input$apply_filters, {
    filtered_data()
  })
  
  observeEvent(input$apply_filters, {
    if (!tabsInserted()) {
      
      # 4A) MLE Phylogenetic Tree
      insertTab(
        "main_tabs",
        tabPanel(
          "MLE Phylogenetic Tree",
          fluidPage(
            plotlyOutput("phylo_tree"),
            br(),
            helpText("Phylogenetic tree for species: ", input$table_species_filter)
          )
        ),
        target   = "Allele Frequency Plot",
        position = "after"
      )
      
      # 4B) Sequence Alignment (Modified)
      insertTab(
        "main_tabs",
        tabPanel(
          "Sequence Alignment",
          sidebarLayout(
            sidebarPanel(
              selectInput("asv_select", "Select ASVs to Display:",
                          choices = {
                            fasta_file <- "/Users/ronaarqureshi/Desktop/asv.fasta"
                            sequences  <- readDNAStringSet(fasta_file, format = "fasta")
                            names(sequences)
                          },
                          multiple = TRUE)
            ),
            mainPanel(
              htmlOutput("sequence_alignment_html"),
              br(),
              reactableOutput("nucleotide_diff_table")
            )
          )
        ),
        target   = "MLE Phylogenetic Tree",
        position = "after"
      )
      
      # 4C) Geospatial Map
      insertTab(
        "main_tabs",
        tabPanel(
          "Geospatial Map",
          sidebarLayout(
            sidebarPanel(
              selectInput("map_asv_select", "Select ASVs to show on Map:",
                          choices  = unique(metadata()$asv_md$ASV_ID),
                          selected = unique(metadata()$asv_md$ASV_ID),
                          multiple = TRUE),
              
              radioButtons("map_filter_mode", "Map Filter Mode:",
                           choices  = c("Show All" = "all",
                                        "Unique ASVs" = "unique",
                                        "Common ASVs" = "common",
                                        "Custom Filter" = "custom"),
                           selected = "all"),
              
              conditionalPanel(
                condition = "input.map_filter_mode == 'custom'",
                selectInput("custom_region_select", "Select Regions:",
                            choices  = c("crete", "thessaly", "western macedonia", "central macedonia"),
                            multiple = TRUE,
                            selected = c("crete", "thessaly")),
                
                radioButtons("custom_filter_type", "Custom Filter Type:",
                             choices  = c("Common", "Unique"),
                             selected = "Common")
              )
            ),
            mainPanel(
              leafletOutput("greece_map", height = "600px"),
              br(),
              reactableOutput("map_asv_table")
            )
          )
        ),
        target   = "Sequence Alignment",
        position = "after"
      )
      
      tabsInserted(TRUE)
    }
  })
  
  # "Return to Filtering" button
  output$return_button_ui <- renderUI({
    if (tabsInserted()) {
      actionButton("return_to_filter", "Return to Filtering")
    }
  })
  observeEvent(input$return_to_filter, {
    removeTab("main_tabs", target = "MLE Phylogenetic Tree")
    removeTab("main_tabs", target = "Sequence Alignment")
    removeTab("main_tabs", target = "Geospatial Map")
    tabsInserted(FALSE)
    updateTabsetPanel(session, "main_tabs", selected = "Table View")
  })
  
  
  output$phylo_tree <- renderPlotly({
    req(new_dataset())
    req(input$table_species_filter)
    
    species <- input$table_species_filter
    asv_ids <- new_dataset()$md$ASV_ID
    
    fasta_file <- "/Users/ronaarqureshi/Desktop/asv.fasta"
    sequences  <- readDNAStringSet(fasta_file, format = "fasta")
    sequences  <- sequences[names(sequences) %in% as.character(asv_ids)]
    
    if (length(sequences) < 2) {
      return(plotly::plotly_empty())
    }
    
    alignment       <- msa(sequences)
    aligned_phang   <- msaConvert(alignment, type = "phangorn")
    dna_bin         <- as.DNAbin(aligned_phang)
    dist_matrix     <- dist.dna(dna_bin, model = "raw")
    tree            <- nj(dist_matrix)
    tree$tip.label  <- names(sequences)
    
    # For demonstration, random region assignment
    tip_ids  <- tree$tip.label
    locations <- c("Crete", "thessaly", "western macedonia", "central macedonia")
    set.seed(42)
    grp <- sample(locations, size = length(tip_ids), replace = TRUE)
    dat <- data.frame(id = tip_ids, grp = grp, stringsAsFactors = FALSE)
    
    p1    <- ggtree(tree)
    metat <- p1$data %>% dplyr::inner_join(dat, by = c("label" = "id"))
    
    p2 <- p1 +
      geom_point(
        data = metat,
        aes(
          x = x, y = y, colour = label,
          text = paste("Region:", grp)
        )
      ) +
      labs(
        title = paste("Interactive Phylogenetic Tree for", species),
        color = "ASV"
      )
    
    ggplotly(p2, tooltip = "text")
  })
  
  # 6. Sequence Alignment (Modified)
  output$sequence_alignment_html <- renderUI({
    fasta_file <- "/Users/ronaarqureshi/Desktop/asv.fasta"
    dna_sequences  <- readDNAStringSet(fasta_file, format = "fasta")
    
    selected_asvs <- input$asv_select
    if (is.null(selected_asvs) || length(selected_asvs) == 0) {
      return(HTML("<pre>No ASVs selected.</pre>"))
    }
    
    # Filter the DNA sequences based on selection
    filtered_dna <- dna_sequences[names(dna_sequences) %in% selected_asvs]
    if (length(filtered_dna) == 0) {
      return(HTML("<pre>No matching ASVs found in the fasta file.</pre>"))
    }
    
    # Translate DNA to amino acid sequences
    aa_sequences <- translate(filtered_dna)
    
    # Align amino acid sequences
    aa_alignment <- msa(aa_sequences)
    aligned_aa <- as.character(aa_alignment)
    # Convert alignment to a matrix (each sequence as a vector of characters)
    aligned_matrix <- do.call(rbind, strsplit(aligned_aa, split = ""))
    rownames(aligned_matrix) <- names(filtered_dna)
    
    # Compute consensus sequence for amino acids
    consensus <- apply(aligned_matrix, 2, function(col) {
      ux <- unique(col)
      ux[which.max(tabulate(match(col, ux)))]
    })
    
    consensus_str <- paste(consensus, collapse = "")
    
    # Build HTML output for alignment with highlighted differences
    alignment_html <- "<pre>"
    alignment_html <- paste0(alignment_html, "Consensus: ", consensus_str, "\n\n")
    
    for(i in 1:nrow(aligned_matrix)) {
      seq_label <- rownames(aligned_matrix)[i]
      alignment_html <- paste0(alignment_html, seq_label, ": ")
      seq_chars <- aligned_matrix[i, ]
      for(j in seq_along(seq_chars)) {
        char <- seq_chars[j]
        if(char != consensus[j]) {
          # Highlight differences in red and label with the position number as a tooltip
          alignment_html <- paste0(alignment_html, "<span style='color:red;' title='Position ", j, "'>", char, "</span>")
        } else {
          alignment_html <- paste0(alignment_html, char)
        }
      }
      alignment_html <- paste0(alignment_html, "\n")
    }
    alignment_html <- paste0(alignment_html, "</pre>")
    HTML(alignment_html)
  })
  
  output$nucleotide_diff_table <- renderReactable({
    fasta_file <- "/Users/ronaarqureshi/Desktop/asv.fasta"
    dna_sequences  <- readDNAStringSet(fasta_file, format = "fasta")
    
    selected_asvs <- input$asv_select
    if (is.null(selected_asvs) || length(selected_asvs) == 0) {
      return(NULL)
    }
    
    # Filter the DNA sequences based on selection
    filtered_dna <- dna_sequences[names(dna_sequences) %in% selected_asvs]
    if (length(filtered_dna) < 2) {
      return(NULL)
    }
    
    # Align nucleotide sequences
    dna_alignment <- msa(filtered_dna)
    aligned_dna <- as.character(dna_alignment)
    dna_matrix <- do.call(rbind, strsplit(aligned_dna, split = ""))
    rownames(dna_matrix) <- names(filtered_dna)
    
    n <- nrow(dna_matrix)
    diff_matrix <- matrix(0, n, n)
    rownames(diff_matrix) <- rownames(dna_matrix)
    colnames(diff_matrix) <- rownames(dna_matrix)
    
    for(i in 1:(n-1)) {
      for(j in (i+1):n) {
        diff_count <- sum(dna_matrix[i, ] != dna_matrix[j, ])
        diff_matrix[i, j] <- diff_count
        diff_matrix[j, i] <- diff_count
      }
    }
    
    # Convert difference matrix to a data frame for display
    diff_df <- as.data.frame(diff_matrix)
    diff_df <- cbind(ASV = rownames(diff_df), diff_df)
    
    reactable(diff_df, columns = list(
      ASV = colDef(name = "ASV")
      # The remaining columns show the pairwise nucleotide differences
    ))
  })
  
  # 7. Geospatial Map
  map_data <- reactive({
    req(metadata()$asv_read, metadata()$sample_md)
    
    cat("\n--- DEBUG: Building map_data() ---\n")
    cat("DEBUG: unique region names in metadata()$sample_md$REGION_L1:\n")
    print(unique(metadata()$sample_md$REGION_L1))
    
    selected_asvs <- input$map_asv_select
    if (is.null(selected_asvs) || length(selected_asvs) == 0) {
      selected_asvs <- filtered_data()$md$ASV_ID
    }
    
    cat("DEBUG: user selected these ASVs for the map:\n")
    print(selected_asvs)
    
    asv_long <- metadata()$asv_read %>%
      pivot_longer(cols = -SAMPLE_ID, names_to = "ASV_ID", values_to = "READ_NUM") %>%
      mutate(ASV_ID = as.character(ASV_ID)) %>%
      filter(ASV_ID %in% selected_asvs, READ_NUM > 0)
    
    cat("DEBUG: after pivot/filter, asv_long has rows:\n")
    print(nrow(asv_long))
    
    asv_long <- asv_long %>%
      left_join(metadata()$sample_md, by = "SAMPLE_ID") %>%
      mutate(region = tolower(trimws(as.character(REGION_L1))))
    
    cat("DEBUG: unique region values in asv_long$region:\n")
    print(unique(asv_long$region))
    
    target_regions <- c("crete", "thessaly", "western macedonia", "central macedonia")
    asv_long <- asv_long %>% filter(region %in% target_regions)
    
    cat("DEBUG: rows left after filtering to target regions:\n")
    print(nrow(asv_long))
    
    mode <- input$map_filter_mode
    if (mode == "all") {
      filtered_asv_long <- asv_long
      custom_regions    <- target_regions
    } else if (mode == "unique") {
      asv_region_counts <- asv_long %>%
        group_by(ASV_ID) %>%
        summarize(reg_count = n_distinct(region), .groups = "drop")
      unique_asvs <- asv_region_counts %>% filter(reg_count == 1) %>% pull(ASV_ID)
      filtered_asv_long <- asv_long %>% filter(ASV_ID %in% unique_asvs)
      custom_regions    <- target_regions
    } else if (mode == "common") {
      asv_region_counts <- asv_long %>%
        group_by(ASV_ID) %>%
        summarize(reg_count = n_distinct(region), .groups = "drop")
      common_asvs <- asv_region_counts %>% filter(reg_count == length(target_regions)) %>% pull(ASV_ID)
      filtered_asv_long <- asv_long %>% filter(ASV_ID %in% common_asvs)
      custom_regions    <- target_regions
    } else if (mode == "custom") {
      custom_regions <- input$custom_region_select
      if (is.null(custom_regions) || length(custom_regions) == 0) {
        custom_regions <- target_regions
      }
      filtered_asv_long <- asv_long %>% filter(region %in% custom_regions)
      
      if (input$custom_filter_type == "Common") {
        asv_region_counts <- filtered_asv_long %>%
          group_by(ASV_ID) %>%
          summarize(reg_count = n_distinct(region), .groups = "drop")
        common_asvs <- asv_region_counts %>% filter(reg_count == length(custom_regions)) %>% pull(ASV_ID)
        filtered_asv_long <- filtered_asv_long %>% filter(ASV_ID %in% common_asvs)
      } else if (input$custom_filter_type == "Unique") {
        asv_region_counts <- filtered_asv_long %>%
          group_by(ASV_ID) %>%
          summarize(reg_count = n_distinct(region), .groups = "drop")
        unique_asvs <- asv_region_counts %>% filter(reg_count == 1) %>% pull(ASV_ID)
        filtered_asv_long <- filtered_asv_long %>% filter(ASV_ID %in% unique_asvs)
      }
    }
    
    cat("DEBUG: final filtered_asv_long rows:\n")
    print(nrow(filtered_asv_long))
    
    region_summary <- filtered_asv_long %>%
      group_by(region) %>%
      summarize(
        count = n_distinct(ASV_ID),
        asvs  = paste(unique(ASV_ID), collapse = ", "),
        .groups = "drop"
      )
    
    region_coords <- data.frame(
      region = c("crete", "thessaly", "western macedonia", "central macedonia"),
      lat    = c(35.2, 39.1, 40.3, 40.7),
      lng    = c(24.8, 22.9, 21.8, 22.9),
      stringsAsFactors = FALSE
    )
    
    region_coords <- region_coords %>%
      filter(region %in% unique(filtered_asv_long$region))
    
    merged_data <- merge(region_coords, region_summary, by = "region", all.x = TRUE)
    merged_data$count[is.na(merged_data$count)] <- 0
    merged_data$asvs[is.na(merged_data$asvs)]   <- "None"
    
    cat("DEBUG: final merged_data for Leaflet:\n")
    print(merged_data)
    
    merged_data
  })
  
  output$greece_map <- renderLeaflet({
    req(map_data())
    data <- map_data()
    
    if (nrow(data) == 0) {
      leaflet() %>%
        addTiles() %>%
        setView(lng = 23.0, lat = 38.0, zoom = 6)
    } else {
      leaflet(data) %>%
        addTiles() %>%
        setView(lng = 23.0, lat = 38.0, zoom = 6) %>%
        addMarkers(~lng, ~lat,
                   popup = ~paste0("<strong>", region, "</strong><br>",
                                   count, " ASVs: ", asvs))
    }
  })
  
  output$map_asv_table <- renderReactable({
    req(map_data())
    reactable(map_data(), columns = list(
      region = colDef(name = "Region"),
      count  = colDef(name = "Number of ASVs"),
      asvs   = colDef(name = "ASVs")
    ))
  })
}

shinyApp(ui, server)
