# Required Libraries
library(readxl)
library(dplyr)
library(igraph)
library(tidyverse)
library(pheatmap)
library(mixOmics)

# Read Metabolomics data
metabo_df <- read.csv("Metabolomics.csv", row.names = 1)
Y <- metabo_df$label
metabo_df <- subset(metabo_df, select = 16:23)

# Read Lipidomics data
lipido_df <- read.csv("Lipidomics.csv", row.names = 1)
lipido_df <- subset(lipido_df, select = 12:19)

# Read Proteomics data
proteo_df <- read.csv("Proteomics.csv", row.names = 1)
proteo_df <- subset(proteo_df, select = 6:13)

# Verify the dimensions
dim(metabo_df)
dim(lipido_df)
dim(proteo_df)

# Verify rownames
rownames(metabo_df)
rownames(lipido_df)
rownames(proteo_df)

# Prepare the data list for DIABLO
X <- list(metabolomics = metabo_df, 
          lipidomics = lipido_df, 
          proteomics = proteo_df)

# Assuming you have a design matrix
design = matrix(0.1, ncol = 3, nrow = 3, 
                dimnames = list(c("metabolomics", "lipidomics", "proteomics"), 
                                c("metabolomics", "lipidomics", "proteomics")))
diag(design) = 0

# Define the keepX list to specify the number of variables to retain
keepX <- list(metabolomics = c(5, 5), lipidomics = c(5, 5), proteomics = c(5, 5))

# Run block.splsda with specified keepX
result.splsda <- block.splsda(X, Y, ncomp = 2,
                              keepX = keepX, design = design)

# Plot DIABLO results with 1 component, without legend
png(filename = "plotDiablo.png", width = 800, height = 600)
plotDiablo(result.splsda, ncomp = 1, legend = TRUE)
dev.off()

# Plot the individual samples without names but with a legend
png(filename = "plotIndiv.png", width = 800, height = 600)
plotIndiv(result.splsda, ind.names = FALSE, legend = TRUE, 
          title = 'Individual plot, DIABLO comp 1 - 2')
dev.off()

# Plot the variables with custom style and legend
png(filename = "plotVar.png", width = 800, height = 600)
plotVar(result.splsda, var.names = FALSE, style = 'graphics', legend = TRUE, 
        pch = c(16, 17, 15), cex = c(2,2,2), 
        col = c('darkorchid', 'brown1', 'lightgreen'),
        title = 'Variables plot, DIABLO comp 1 - 2')
dev.off()

# Plot arrows to show the relationship between components and classes
png(filename = "plotArrow.png", width = 800, height = 600)
plotArrow(result.splsda, ind.names = FALSE, legend = TRUE, 
          title = 'Arrow plot, DIABLO comp 1 - 2')
dev.off()

# Create a circos plot with specified cutoff and colors
png(filename = "circosPlot.png", width = 800, height = 800)
circosPlot(result.splsda, cutoff = 0.7, line = TRUE, 
           color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
           color.cor = c("chocolate3", "grey20"), size.labels = 1.5)
dev.off()

# Create a network plot
network(result.splsda, blocks = c(1, 2, 3), 
        cutoff = 0.4, size.node = 0.08,
        color.node = c('darkorchid', 'brown1', 'lightgreen'),
        save = 'png', name.save = 'diablo-network-plot')

# Plot loadings for the first component
png(filename = "plotLoadings.png", width = 800, height = 600)
plotLoadings(result.splsda, comp = 1, contrib = 'max', method = 'median')
dev.off()

# Create CIM plot
png(filename = "cimDiablo_plot.png", width = 800, height = 800)
cimDiablo(result.splsda, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = 1, margin=c(8, 20), legend.position = "topright", trim = TRUE)
dev.off()
