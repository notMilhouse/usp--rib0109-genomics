library(ggplot2)
library(dplyr)
library(pheatmap)

# 1. Distribution of Gene Expression (Density Plot)
ggplot(rnaseq_data, aes(x = log10(jovem_count + 1), color = "Jovem")) +
  geom_density() +
  geom_density(aes(x = log10(naojovem_count + 1), color = "NaoJovem")) +
  scale_color_manual(values = c("blue", "red"), name = "Category") +
  labs(title = "Gene Expression Distribution (Log Scale)", x = "Log10(Gene Expression + 1)", y = "Density") +
  theme_minimal()

# 2. Boxplot of Gene Expression Between Groups
ggplot(rnaseq_data, aes(x = "Jovem", y = log10(jovem_count + 1))) +
  geom_boxplot(fill = "blue", alpha = 0.6) +
  geom_boxplot(aes(x = "NaoJovem", y = log10(naojovem_count + 1)), fill = "red", alpha = 0.6) +
  labs(title = "Gene Expression Boxplot", x = "Group", y = "Log10(Gene Expression + 1)") +
  theme_minimal()

# 3. Top Expressed Genes Bar Plot
top_genes <- rnaseq_data %>%
  mutate(total_count = jovem_count + naojovem_count) %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 10)

ggplot(top_genes, aes(x = reorder(gene_id, -total_count), y = total_count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Top 10 Expressed Genes", x = "Gene ID", y = "Total Count") +
  theme_minimal() +
  coord_flip()

# 4. Correlation Analysis
correlation <- cor(rnaseq_data$jovem_count, rnaseq_data$naojovem_count, method = "pearson")
cat("Pearson Correlation between Jovem and NaoJovem gene expression:", correlation, "\n")

ggplot(rnaseq_data, aes(x = jovem_count, y = naojovem_count)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(title = "Correlation Between Jovem and NaoJovem Expression", x = "Jovem Count", y = "NaoJovem Count") +
  theme_minimal()

# 5. Heatmap of Top Genes Expression
# Select top 20 genes for heatmap
heatmap_genes <- rnaseq_data %>%
  mutate(total_count = jovem_count + naojovem_count) %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 20) %>%
  select(gene_id, jovem_count, naojovem_count)

# Prepare data for heatmap
heatmap_matrix <- as.matrix(heatmap_genes[, -1])  # Exclude gene_id
rownames(heatmap_matrix) <- heatmap_genes$gene_id

# Generate heatmap
pheatmap(
  log10(heatmap_matrix + 1), 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  color = colorRampPalette(c("navy", "white", "firebrick"))(50), 
  main = "Heatmap of Top 20 Expressed Genes"
)

# Identify top 2 genes based on total expression
rnaseq_data$total_expression <- rnaseq_data$jovem_count + rnaseq_data$naojovem_count
top_genes <- rnaseq_data[order(-rnaseq_data$total_expression), ][1:2, ]

# Extract top genes
top_gene_1 <- top_genes[1, ]
top_gene_2 <- top_genes[2, ]

# Perform Wilcoxon tests for differential expression
gene_1_test <- wilcox.test(c(top_gene_1$jovem_count), c(top_gene_1$naojovem_count))
gene_2_test <- wilcox.test(c(top_gene_2$jovem_count), c(top_gene_2$naojovem_count))

# Print test results
cat("Top Gene 1:", top_gene_1$gene_id, "\n")
cat("Wilcoxon p-value:", gene_1_test$p.value, "\n\n")
cat("Top Gene 2:", top_gene_2$gene_id, "\n")
cat("Wilcoxon p-value:", gene_2_test$p.value, "\n\n")

# Prepare data for boxplots
top_gene_1_long <- data.frame(
  Category = c(rep("jovem", length(top_gene_1$jovem_count)),
               rep("naojovem", length(top_gene_1$naojovem_count))),
  Expression = c(top_gene_1$jovem_count, top_gene_1$naojovem_count)
)

top_gene_2_long <- data.frame(
  Category = c(rep("jovem", length(top_gene_2$jovem_count)),
               rep("naojovem", length(top_gene_2$naojovem_count))),
  Expression = c(top_gene_2$jovem_count, top_gene_2$naojovem_count)
)

# Create boxplots
plot_1 <- ggplot(top_gene_1_long, aes(x = Category, y = Expression, fill = Category)) +
  geom_boxplot() +
  labs(title = paste("Expression of", top_gene_1$gene_id),
       x = "Category", y = "Expression Count") +
  theme_minimal()

plot_2 <- ggplot(top_gene_2_long, aes(x = Category, y = Expression, fill = Category)) +
  geom_boxplot() +
  labs(title = paste("Expression of", top_gene_2$gene_id),
       x = "Category", y = "Expression Count") +
  theme_minimal()

# Display plots
print(plot_1)
print(plot_2)


# Calculate mean expression difference for each gene
rnaseq_data$difference <- abs(rnaseq_data$jovem_count - rnaseq_data$naojovem_count)

# Sort by the biggest difference
top_genes <- rnaseq_data[order(-rnaseq_data$difference), ]
top_genes <- head(top_genes, 5)  # Select top 5 genes with the largest differences

# Prepare data for plotting
long_data <- data.frame(
  gene_id = rep(top_genes$gene_id, each = 2),
  Category = rep(c("jovem", "naojovem"), times = nrow(top_genes)),
  Expression = c(top_genes$jovem_count, top_genes$naojovem_count)
)

# Plot the top genes
plot <- ggplot(long_data, aes(x = gene_id, y = Expression, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Genes with Largest Expression Differences",
       x = "Gene", y = "Expression Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot)

# Optional: Save the plot
ggsave("largest_expression_diff_genes.png", plot, width = 10, height = 6)

# Print table for reference
top_genes
