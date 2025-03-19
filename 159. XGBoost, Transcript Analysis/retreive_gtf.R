# import library
library(GenomicFeatures)

# Function to download the gtf file
download_gtf <- function(url, filename) {
  if (!file.exists(filename)) {
    download.file(url, destfile = filename, method = "wget")
  } 
  else {
    message("File '", filename, "' already exists. Skipping download.")
  }
}

# Download the gtf file
gtf_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz"
gtf_filename <- "gencode.v45.primary_assembly.annotation.gtf.gz"
download_gtf(gtf_url, gtf_filename)

# Create TxDb object from gtf file
txdb <- makeTxDbFromGFF(gtf_filename, format = "gtf")

# Build S4 object containing all genes and transcript IDs
gene_transcripts <- select(txdb, keys(txdb, keytype = "GENEID"), "TXNAME", "GENEID")

# Compute mean, minimum, and maximum number of transcripts per gene
num_transcripts <- table(gene_transcripts$GENEID)
mean_num_transcripts <- mean(num_transcripts)
min_num_transcripts <- min(num_transcripts)
max_num_transcripts <- max(num_transcripts)

# Plot the histogram as a bar plot without gene names on the x-axis
barplot(num_transcripts, 
        main = "Histogram of Transcript Counts per Gene",
        xlab = "Genes",
        ylab = "Number of Transcripts",
        col = "skyblue",
        border = "black",
        names.arg = rep("", length(num_transcripts)))

# Save S4 object as .rds file
saveRDS(gene_transcripts, file = "gene_transcripts.rds")
