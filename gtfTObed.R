# Read the GTF-like data
gtf_data <- read.table("Sorghum_bicolor.Sorghum_bicolor_NCBIv3.57.chr.gtf", sep = "\t", header = FALSE, stringsAsFactors = FALSE)

# Filter to include only 'exon' and 'CDS' features
filtered_data <- subset(gtf_data, gtf_data$V3 == "exon" | gtf_data$V3 == "CDS")

# Extract gene IDs from attributes
gene_ids <- character(nrow(filtered_data))

for (i in 1:nrow(filtered_data)) {
  attributes <- strsplit(as.character(filtered_data[i, 9]), ";")[[1]]
  gene_id <- grep("gene_id", attributes, value = TRUE)
  
  if (length(gene_id) > 0) {
    gene_ids[i] <- sub(".*gene_id \"([^\"]+)\".*", "\\1", gene_id)
  } else {
    gene_ids[i] <- NA
  }
}

# Create a data frame with BED-like columns and include gene IDs
genes_bed <- data.frame(
  chrom = filtered_data$V1,       # Chromosome column
  chromStart = filtered_data$V4,  # Start position column
  chromEnd = filtered_data$V5,    # End position column
  strand = ".",                   # Strand column (if not available)
  score = ".",                    # Score column (if not available)
  gene_id = gene_ids              # Gene IDs extracted from attributes
)

# Write the data to a BED file
write.table(genes_bed, file = "genes_with_ids.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
