setwd("C:\\Users\\admin\\College\\College\\Coding\\gse19188_analyze")
library(affy)
library(GEOquery)
library(tidyverse)
library(dplyr)
getGEOSuppFiles("GSE19188")
untar("GSE19188/GSE19188_RAW.tar", exdir = "data/")
raw.data <- ReadAffy(celfile.path = "data/")
normalized.data <- rma(raw.data)
normalized.expr <- as.data.frame(exprs(normalized.data))
gse <- getGEO("GSE19188", GSEMatrix = TRUE)
feature.data <- gse$GSE19188_series_matrix.txt.gz@featureData@data
feature.gene$ID <- feature.data$ID
feature.gene$`Gene Symbol` <- feature.data$`Gene Symbol`
feature.gene <- filter(feature.gene, !is.na("Gene Symbol"))
normalized.expr <- normalized.expr %>%
  rownames_to_column(var = 'ID') %>%
  left_join(
    feature.gene, 
    by = 'ID'
  )

normalized.expr <- normalized.expr %>%
  select(ID,`Gene Symbol`, everything())
normalized.expr <- filter(normalized.expr, `Gene Symbol` != "")
normalized.expr$mean_val <- rowMeans(select(normalized.expr, where(is.numeric)))
clean.data <- normalized.expr %>%
  group_by(`Gene Symbol`) %>%
  slice_max(order_by = mean_val, n = 1, with_ties = FALSE) %>%
  ungroup()

pheno.data <- gse$GSE19188_series_matrix.txt.gz@phenoData@data
metadata <- pheno.data %>%
  select(geo_accession, characteristics_ch1, characteristics_ch1.1, characteristics_ch1.2, characteristics_ch1.3, characteristics_ch1.4, title)
metadata <- metadata %>%
  mutate(
    condition = case_when(
      grepl("healthy", characteristics_ch1, ignore.case = TRUE) ~ "healthy",
      grepl("tumor", characteristics_ch1, ignore.case = TRUE) ~ "cancer",
      TRUE ~ "unknown"
    )
  )
metadata$characteristics_ch1.1 <- gsub("cell type: ", "", metadata$characteristics_ch1.1)
metadata$characteristics_ch1 <- NULL
metadata <- metadata %>%
  select(geo_accession,
         cell_type = characteristics_ch1.1,
         survival = characteristics_ch1.2 , 
         survival_status = characteristics_ch1.3,
         gender = characteristics_ch1.4,
         sample_name = title,
         condition)
metadata$survival <- gsub("overall survival: ", "", metadata$survival)
metadata$survival_status <- gsub("status: ", "", metadata$survival_status)
metadata$gender <- gsub("gender: ", "", metadata$gender)
dir.create("processed_data")
write.csv(metadata,
          "processed_data/metadata.csv",
          row.names = FALSE)
#normalized.expr
write.csv(clean.data,
          "processed_data/normalized_expression.csv",
          row.names = FALSE)