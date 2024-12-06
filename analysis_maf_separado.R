library(ggplot2)
library(dplyr)

# Separate mutation frequencies for jovem and naojovem
mutation_counts_separated <- maf_data %>%
  group_by(Hugo_Symbol, Category) %>%
  summarize(Frequency = n(), .groups = 'drop')

# Get top 10 mutations for each category
top_mutations_jovem <- mutation_counts_separated %>%
  filter(Category == "jovem") %>%
  arrange(desc(Frequency)) %>%
  slice_head(n = 10)

top_mutations_naojovem <- mutation_counts_separated %>%
  filter(Category == "naojovem") %>%
  arrange(desc(Frequency)) %>%
  slice_head(n = 10)

# Combine top mutations for plotting
combined_top_mutations <- mutation_counts_separated %>%
  filter(Hugo_Symbol %in% c(top_mutations_jovem$Hugo_Symbol, top_mutations_naojovem$Hugo_Symbol))

# Plot jovem mutation frequencies
plot_mutation_freq_jovem <- ggplot(combined_top_mutations %>% filter(Category == "jovem"),
                                   aes(x = Hugo_Symbol, y = Frequency)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Top 10 Mutated Genes (Jovem)",
       x = "Gene", y = "Mutation Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_mutation_freq_jovem)

# Plot naojovem mutation frequencies
plot_mutation_freq_naojovem <- ggplot(combined_top_mutations %>% filter(Category == "naojovem"),
                                      aes(x = Hugo_Symbol, y = Frequency)) +
  geom_bar(stat = "identity", fill = "salmon") +
  labs(title = "Top 10 Mutated Genes (Naojovem)",
       x = "Gene", y = "Mutation Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_mutation_freq_naojovem)

# Save plots
ggsave("mutation_frequencies_jovem.png", plot_mutation_freq_jovem, width = 10, height = 6)
ggsave("mutation_frequencies_naojovem.png", plot_mutation_freq_naojovem, width = 10, height = 6)

# Identify overlapping genes between the two groups
overlapping_genes <- intersect(top_mutations_jovem$Hugo_Symbol, top_mutations_naojovem$Hugo_Symbol)
cat("Overlapping genes between Jovem and Naojovem:\n")
print(overlapping_genes)
