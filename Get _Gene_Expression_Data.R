library(affy)
library(GEOquery)
library(limma)
library(tidyverse)

gene_expression_data <- ReadAffy(celfile.path = "D:/UTM/20222023 SEM 2/PSM I/FYP/dataset/Gene_Expression")

gene_expression_data <- rma(gene_expression_data)

gene_expression_data <- as.data.frame(exprs(gene_expression_data))

# normalizing and calculating expression
gse <- getGEO("GSE20347", GSEMatrix = TRUE)

# join the variable and sample into table
feature_data <- gse$GSE20347_series_matrix.txt.gz@featureData
feature_data <- feature_data[ , c(1, 11)]

gene_expression_data <-gene_expression_data %>%
  rownames_to_column(var='ID')%>%
  inner_join(., feature_data@data, by = 'ID')

# rearrange the table
gene_expression_data <- gene_expression_data %>%
  dplyr::select("Gene Symbol", c(2:35))

write.csv(gene_expression_data, "D:/UTM/20222023 SEM 2/PSM I/FYP/dataset/GeneExpression.csv", row.names = FALSE)
