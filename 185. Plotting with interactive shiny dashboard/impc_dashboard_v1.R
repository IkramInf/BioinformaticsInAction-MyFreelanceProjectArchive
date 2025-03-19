# Load the Required Libraries
library(shiny)
library(shinydashboard)
library(tidyverse)
library(plotly)
library(cluster)

ui <- dashboardPage(
  dashboardHeader(title = "Mouse Phenotypic Analysis Dashboard"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Phenotype Analysis", tabName = "phenotype"),
      menuItem("Genotype Analysis", tabName = "genotype"),
      menuItem("Gene Clusters", tabName = "clusters")),
    selectInput("phenotype_var", "Select Phenotype:", choices = NULL),
    numericInput("top_n", "Top N genes:", value = 10, min = 5, max = 100),
    numericInput("n_clusters", "Number of Clusters:", value = 3, min = 2, max = 10)),

  dashboardBody(
    tabItems(
      tabItem(tabName = "phenotype",
              fluidRow(box(plotlyOutput("phenotype_plot"), width = 12))),
      
      tabItem(tabName = "genotype",
              fluidRow(box(plotlyOutput("genotype_plot"), width = 12))),
      
      tabItem(tabName = "clusters",
              fluidRow(box(plotlyOutput("cluster_plot"), width = 12, height = 600)))
    )
  )
)

server <- function(input, output, session) {
  # Read and arrange data
  data <- reactive({read.csv("cleaned_data.csv")})
  observe({updateSelectInput(session, "phenotype_var",
                      choices = unique(data()$parameter_name))})
  
  # Output phenotype plot
  output$phenotype_plot <- renderPlotly({
    plot_data <- data() %>%
      group_by(parameter_name) %>%
      summarise(mean_pvalue = mean(pvalue),
          n_genes = n_distinct(gene_symbol)) %>%
      arrange(mean_pvalue) %>%
      head(input$top_n)
    
    p <- ggplot(plot_data, 
                aes(x = reorder(parameter_name, -mean_pvalue),
                    y = mean_pvalue,
                    text = paste("Genes:", n_genes))) +
      geom_bar(stat = "identity", fill = "lightblue") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Phenotype", y = "Mean p-value")
    
    ggplotly(p, tooltip = c("x", "y", "text"))
  })
  
  # Output genotype plot
  output$genotype_plot <- renderPlotly({
    req(input$phenotype_var)
    
    plot_data <- data() %>%
      filter(parameter_name == input$phenotype_var) %>%
      group_by(gene_symbol) %>%
      summarise(
        mean_pvalue = mean(pvalue, na.rm = TRUE),
        n_measurements = n()
      ) %>%
      arrange(mean_pvalue) %>%
      head(input$top_n)
    
    p <- ggplot(plot_data, 
                aes(x = reorder(gene_symbol, mean_pvalue),
                    y = mean_pvalue,
                    text = paste("Gene:", gene_symbol,
                                 "\nP-value:", round(mean_pvalue, 4),
                                 "\nMeasurements:", n_measurements))) +
      geom_point(size = 3, color = "darkred") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(x = "Gene", 
           y = "P-value",
           title = paste("Top Genes for", input$phenotype_var)) +
      coord_flip()
    
    ggplotly(p, tooltip = "text")
  })
  
  # Output gene clusters plot
  clusters <- reactive({
    matrix_data <- data() %>%
      group_by(gene_symbol, parameter_name) %>%
      summarise(mean_pvalue = mean(pvalue), .groups = 'drop') %>%
      pivot_wider(names_from = parameter_name,
                  values_from = mean_pvalue,
                  values_fill = 1)
    
    genes <- matrix_data$gene_symbol
    matrix_data <- as.matrix(matrix_data[,-1])
    
    pca <- prcomp(matrix_data, scale. = TRUE)
    coords <- pca$x[,1:2]
    
    km <- kmeans(matrix_data, centers = input$n_clusters)
    
    tibble(gene = genes, cluster = factor(km$cluster), PC1 = coords[,1], PC2 = coords[,2])
  })
  
  output$cluster_plot <- renderPlotly({
    p <- ggplot(clusters(), aes(x = PC1, y = PC2, color = cluster, text = gene)) +
      geom_point(size = 3, alpha = 0.6) +
      theme_minimal() +
      labs(x = "PC1", y = "PC2", title = "Gene Clusters")
    
    ggplotly(p, tooltip = c("text", "color"))
  })
}

# Run the App
shinyApp(ui = ui, server = server)
