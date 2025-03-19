# Load Required Libraries
library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(cluster)
library(DBI)
library(RMySQL)

# Database connection function
establish_db_connection <- function() {
  dbConnect(RMySQL::MySQL(),
            dbname = "impcdb1",
            host = "localhost",
            user = "root",
            password = ""
  )
}

ui <- dashboardPage(
  dashboardHeader(title = "Mouse Phenotypic Analysis Dashboard"),
  
  dashboardSidebar(
    selectInput("analysis_option", "Select Analysis:",
                choices = c("SelectPhenotypeViewGenes", "SelectGeneViewPhenotypes",
                            "MouseLifeStage", "DiseaseClusters")),
    
    # Dynamic inputs based on analysis type
    uiOutput("dynamic_inputs")
  ),
  
  dashboardBody(
    plotlyOutput("main_plot")
  )
)

server <- function(input, output, session) {
  # Create reactive database connection
  con <- reactive({
    establish_db_connection()
  })
  
  # Clean up connection when session ends
  onStop(function() {
    if(!is.null(con())) dbDisconnect(con())
  })
  
  # Reactive queries for data
  mouse_data <- reactive({
    req(con())
    dbGetQuery(con(), "SELECT parameter_name, pvalue, gene_symbol, mouse_life_stage 
               FROM Mouse")
  })
  
  disease_data <- reactive({
    req(con())
    dbGetQuery(con(), "SELECT disease_id, disease_term, phenodigm_score 
               FROM Human_Disease")
  })
  
  # Dynamic UI elements
  output$dynamic_inputs <- renderUI({
    switch(input$analysis_option,
           "SelectPhenotypeViewGenes" = list(
             selectInput("phenotype", "Select Phenotype:", 
                         choices = unique(mouse_data()$parameter_name)),
             numericInput("top_n", "Top N:", value = 10, min = 5, max = 100)),
           "SelectGeneViewPhenotypes" = list(
             selectInput("gene", "Select Gene:", 
                         choices = unique(mouse_data()$gene_symbol)),
             numericInput("top_n", "Top N:", value = 10, min = 5, max = 100)),
           "MouseLifeStage" = selectInput("life_stage", "Select MouseLifeStage:", 
                                          choices = unique(mouse_data()$mouse_life_stage)),
           "DiseaseClusters" = numericInput("n_clusters", "Number of Clusters:", 
                                            value = 3, min = 2, max = 10))
  })
  
  # Main plot body
  output$main_plot <- renderPlotly({
    req(input$analysis_option)
    
    plot_data <- switch(input$analysis_option,
                        "SelectPhenotypeViewGenes" = {
                          req(input$phenotype, input$top_n)
                          # First get top N genes by any p-value
                          top_genes <- mouse_data() %>%
                            filter(parameter_name == input$phenotype) %>%
                            slice_min(order_by = pvalue, n = input$top_n) %>%
                            pull(gene_symbol) %>%
                            unique()
                          
                          # Then get ALL p-values for these genes
                          mouse_data() %>%
                            filter(parameter_name == input$phenotype,
                                   gene_symbol %in% top_genes) %>%
                            mutate(significant = pvalue < 0.05)
                        },
                        "SelectGeneViewPhenotypes" = {
                          req(input$gene, input$top_n)
                          # First get top N phenotypes by any p-value
                          top_phenotypes <- mouse_data() %>%
                            filter(gene_symbol == input$gene) %>%
                            slice_min(order_by = pvalue, n = input$top_n) %>%
                            pull(parameter_name) %>%
                            unique()
                          
                          # Then get ALL p-values for these phenotypes
                          mouse_data() %>%
                            filter(gene_symbol == input$gene,
                                   parameter_name %in% top_phenotypes) %>%
                            mutate(significant = pvalue < 0.05)
                        },
                        "MouseLifeStage" = {
                          req(input$life_stage)
                          mouse_data() %>%
                            filter(mouse_life_stage == input$life_stage) %>%
                            group_by(parameter_name) %>%
                            summarise(value = mean(pvalue, na.rm = TRUE))},
                        "DiseaseClusters" = {
                          req(input$n_clusters)
                          disease_clean <- disease_data() %>%
                            group_by(disease_id, disease_term) %>%
                            summarise(value = mean(phenodigm_score), .groups = 'drop')
                          
                          # Create distance matrix and perform PCA
                          dist_matrix <- dist(scale(disease_clean$value))
                          pca_result <- prcomp(as.matrix(dist_matrix), scale. = TRUE)
                          
                          # Perform clustering
                          km <- kmeans(scale(disease_clean$value), centers = input$n_clusters)
                          
                          # Combine results
                          tibble(
                            disease_id = disease_clean$disease_id,
                            disease_term = disease_clean$disease_term,
                            value = disease_clean$value,
                            cluster = as.factor(km$cluster),
                            PC1 = pca_result$x[,1],
                            PC2 = pca_result$x[,2])}
    )
    
    # Create plot based on analysis type
    p <- if(input$analysis_option == "DiseaseClusters") {
      ggplot(plot_data, aes(x = PC1, y = PC2, color = cluster,
                            text = paste("Disease:", disease_term,
                                         "\nPhenodigm Score:", round(value, 2),
                                         "\nCluster:", cluster))) +
        geom_point(size = 3, alpha = 0.7) +
        scale_color_brewer(palette = "Set1") +
        theme_minimal() +
        labs(x = "PC1", y = "PC2", color = "Cluster")
    } else if(input$analysis_option %in% c("SelectPhenotypeViewGenes", "SelectGeneViewPhenotypes")) {
      is_gene_view <- input$analysis_option == "SelectPhenotypeViewGenes"
      id_col <- if(is_gene_view) quo(gene_symbol) else quo(parameter_name)
      
      ggplot(plot_data, 
             aes(x = !!id_col,
                 y = pvalue,
                 group = factor(pvalue),  # Each p-value gets its own bar
                 text = paste(if(is_gene_view) "Gene:" else "Phenotype:", 
                              !!id_col,
                              "\nP-value:", sprintf("%.4e", pvalue),
                              "\nSignificant:", significant))) +
        geom_col(aes(fill = significant),
                 position = position_dodge(width = 0.8),
                 width = 0.7) +
        scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "maroon")) +
        geom_text(aes(label = ifelse(significant, "*", ""),
                      y = pvalue + max(pvalue) * 0.05),
                  position = position_dodge(width = 0.8),
                  size = 5) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = if(is_gene_view) "Gene" else "Phenotype",
             y = "P-value") +
        guides(fill = "none")
    } else {
      ggplot(plot_data, aes(x = parameter_name,
                            y = value,
                            text = paste("Phenotype:", parameter_name,
                                         "\nP-value:", round(value, 4)))) +
        geom_col(fill = "steelblue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Phenotype", y = "P-value")
    }
    
    ggplotly(p, tooltip = "text")
  })
}

shinyApp(ui = ui, server = server)
