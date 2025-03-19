# Load Required Libraries
library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(cluster)

# Specify Input Filepaths
cleaned_data_path <- "final_cleaned_data.csv"
disease_info_path <- "Disease_information.csv"

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
  # Load data
  data <- reactive(read.csv(cleaned_data_path))
  disease_data <- reactive(read.csv(disease_info_path))
  
  # Dynamic UI elements
  output$dynamic_inputs <- renderUI({
    switch(input$analysis_option,
           "SelectPhenotypeViewGenes" = list(
             selectInput("phenotype", "Select Phenotype:", choices = unique(data()$parameter_name)),
             numericInput("top_n", "Top N:", value = 10, min = 5, max = 100)),
           "SelectGeneViewPhenotypes" = list(
             selectInput("gene", "Select Gene:", choices = unique(data()$gene_symbol)),
             numericInput("top_n", "Top N:", value = 10, min = 5, max = 100)),
           "MouseLifeStage" = selectInput("life_stage", "Select MouseLifeStage:", 
                                      choices = unique(data()$mouse_life_stage)),
           "DiseaseClusters" = numericInput("n_clusters", "Number of Clusters:", 
                                             value = 3, min = 2, max = 10))
  })
  
  # Main plot body
  output$main_plot <- renderPlotly({
    req(input$analysis_option)
    
    plot_data <- switch(input$analysis_option,
                        "SelectPhenotypeViewGenes" = {
                          req(input$phenotype, input$top_n)
                          data() %>%
                            filter(parameter_name == input$phenotype) %>%
                            group_by(gene_symbol) %>%
                            summarise(value = mean(pvalue, na.rm = TRUE)) %>%
                            arrange(value) %>%
                            head(input$top_n)},
                        "SelectGeneViewPhenotypes" = {
                          req(input$gene, input$top_n)
                          data() %>%
                            filter(gene_symbol == input$gene) %>%
                            group_by(parameter_name) %>%
                            summarise(value = mean(pvalue, na.rm = TRUE)) %>%
                            arrange(value) %>%
                            head(input$top_n)},
                        "MouseLifeStage" = {
                          req(input$life_stage)
                          data() %>%
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
        labs(x = "PC1", y = "PC2", color = "Cluster")}
    else {
      ggplot(plot_data, aes(x = if(input$analysis_option == "SelectPhenotypeViewGenes") 
        gene_symbol else parameter_name,
        y = value,
        text = paste(if(input$analysis_option == "SelectPhenotypeViewGenes") 
          "Gene:" else "Phenotype:", 
          if(input$analysis_option == "SelectPhenotypeViewGenes") 
            gene_symbol else parameter_name,
          "\nP-value:", round(value, 4)))) +
        geom_col(fill = "steelblue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = if(input$analysis_option == "SelectPhenotypeViewGenes") "Gene" else "Phenotype",
             y = "P-value")}
    
    ggplotly(p, tooltip = "text")
  })
}

shinyApp(ui = ui, server = server)
