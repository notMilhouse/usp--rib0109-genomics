library(ggplot2)
library(dplyr)

# Mutation frequency per gene
mutation_counts <- maf_data %>%
  group_by(Hugo_Symbol, Category) %>%
  summarize(Frequency = n(), .groups = 'drop')

# Top 10 most mutated genes for visualization
top_mutated_genes <- mutation_counts %>%
  group_by(Hugo_Symbol) %>%
  summarize(Total = sum(Frequency), .groups = 'drop') %>%
  arrange(desc(Total)) %>%
  slice_head(n = 10)

# Filter for top genes for plotting
mutation_counts_filtered <- mutation_counts %>%
  filter(Hugo_Symbol %in% top_mutated_genes$Hugo_Symbol)

# Plot mutation frequencies
plot_mutation_freq <- ggplot(mutation_counts_filtered, aes(x = Hugo_Symbol, y = Frequency, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Top 10 Most Mutated Genes by Category",
       x = "Gene", y = "Mutation Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_mutation_freq)

# Variant classification distribution
variant_distribution <- maf_data %>%
  group_by(Variant_Classification, Category) %>%
  summarize(Frequency = n(), .groups = 'drop')

# Plot variant classification distribution
plot_variant_class <- ggplot(variant_distribution, aes(x = Variant_Classification, y = Frequency, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Distribution of Variant Classifications",
       x = "Variant Classification", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_variant_class)

# Save plots
ggsave("mutation_frequencies.png", plot_mutation_freq, width = 10, height = 6)
ggsave("variant_classifications.png", plot_variant_class, width = 10, height = 6)
