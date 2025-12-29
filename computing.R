setwd("C:\\Users\\admin\\College\\College\\Coding\\gse19188_analyze")
expr <- read.csv("processed_data/normalized_expression.csv")
metadata <- read.csv("processed_data/metadata.csv")
normal <- metadata$condition == "healthy"

mean_normal <- rowMeans(expr[3:ncol(expr):158][,normal])
mean_cancer <- rowMeans(expr[3:ncol(expr):158][,!normal])
log2fc <- mean_cancer - mean_normal
expr <- cbind(expr, log2fc)
expr <- expr %>%
  select(ID, gene_symbol, log2fc, mean_val, everything())
p_value <- apply(expr[, 4:ncol(expr)], 1, function (x) {t.test(x[!normal], x[normal])$p.value})
expr <- cbind(expr, p_value)
expr <- expr %>%
  select(ID, gene_symbol, p_value, log2fc, mean_val, everything())
adj_p_value <- p.adjust(p_value, method = "fdr")
expr <- cbind(expr, adj_p_value)
expr$"-log10(padj)" <- -log10(adj_p_value)
expr <- expr %>%
  select(ID, gene_symbol, "-log10(padj)", adj_p_value, p_value, log2fc, mean_val, everything())
write.csv(expr,
          "processed_data/expr_stats.csv",
          row.names = FALSE)