library(ggplot2)
library(dplyr)

# Calculate mutation frequencies for jovem and naojovem
mutation_counts_separated <- maf_data %>%
  group_by(Hugo_Symbol, Category) %>%
  summarize(Frequency = n(), .groups = 'drop')

# --- Least Expressed Genes Analysis ---

# Get least expressed genes for each category
least_mutations_jovem <- mutation_counts_separated %>%
  filter(Category == "jovem") %>%
  arrange(Frequency) %>%
  slice_head(n = 10)

least_mutations_naojovem <- mutation_counts_separated %>%
  filter(Category == "naojovem") %>%
  arrange(Frequency) %>%
  slice_head(n = 10)

# Combine least expressed genes for plotting
combined_least_mutations <- mutation_counts_separated %>%
  filter(Hugo_Symbol %in% c(least_mutations_jovem$Hugo_Symbol, least_mutations_naojovem$Hugo_Symbol))

# Plot least expressed genes for jovem
plot_least_jovem <- ggplot(combined_least_mutations %>% filter(Category == "jovem"),
                           aes(x = Hugo_Symbol, y = Frequency)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Least Mutated Genes (Jovem)",
       x = "Gene", y = "Mutation Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_least_jovem)

# Plot least expressed genes for naojovem
plot_least_naojovem <- ggplot(combined_least_mutations %>% filter(Category == "naojovem"),
                              aes(x = Hugo_Symbol, y = Frequency)) +
  geom_bar(stat = "identity", fill = "salmon") +
  labs(title = "Least Mutated Genes (Naojovem)",
       x = "Gene", y = "Mutation Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(plot_least_naojovem)

# Save plots
ggsave("least_mutations_jovem.png", plot_least_jovem, width = 10, height = 6)
ggsave("least_mutations_naojovem.png", plot_least_naojovem, width = 10, height = 6)

library(ggplot2)
library(dplyr)

# --- Mutation Frequency Calculation ---
mutation_counts_separated <- maf_data %>%
  group_by(Hugo_Symbol, Category) %>%
  summarize(Frequency = n(), .groups = 'drop')

# --- Outlier Detection ---
# Combine jovem and naojovem frequencies for outlier analysis
mutation_counts_combined <- maf_data %>%
  group_by(Hugo_Symbol) %>%
  summarize(Frequency = n(), .groups = 'drop')

# Boxplot to visualize outliers for both jovem and naojovem
outlier_plot_jovem <- ggplot(mutation_counts_separated %>% filter(Category == "jovem"), aes(y = Frequency, x = Category)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(aes(color = Hugo_Symbol), width = 0.2, alpha = 0.6) +
  labs(title = "Mutation Frequencies - Jovem Outliers", y = "Mutation Count", x = "Jovem") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

outlier_plot_naojovem <- ggplot(mutation_counts_separated %>% filter(Category == "naojovem"), aes(y = Frequency, x = Category)) +
  geom_boxplot(fill = "lightcoral") +
  geom_jitter(aes(color = Hugo_Symbol), width = 0.2, alpha = 0.6) +
  labs(title = "Mutation Frequencies - Naojovem Outliers", y = "Mutation Count", x = "Naojovem") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

# Print the plots
print(outlier_plot_jovem)
print(outlier_plot_naojovem)

# Save the plots
ggsave("outlier_mutations_jovem.png", outlier_plot_jovem, width = 10, height = 6)
ggsave("outlier_mutations_naojovem.png", outlier_plot_naojovem, width = 10, height = 6)

# --- Identifying Outlier Genes ---
# Compute the quantiles for outlier detection
q1 <- quantile(mutation_counts_combined$Frequency, 0.25)
q3 <- quantile(mutation_counts_combined$Frequency, 0.75)
iqr <- q3 - q1

lower_bound <- q1 - 1.5 * iqr
upper_bound <- q3 + 1.5 * iqr

# Outliers for both jovem and naojovem
outlier_genes_jovem <- mutation_counts_separated %>%
  filter(Category == "jovem" & (Frequency < lower_bound | Frequency > upper_bound))

outlier_genes_naojovem <- mutation_counts_separated %>%
  filter(Category == "naojovem" & (Frequency < lower_bound | Frequency > upper_bound))

cat("Outlier Genes in Jovem Category:\n")
print(outlier_genes_jovem)

cat("Outlier Genes in Naojovem Category:\n")
print(outlier_genes_naojovem)

# Combine outlier data for plotting and save it
outlier_genes_combined <- bind_rows(outlier_genes_jovem, outlier_genes_naojovem)

# Plot the outlier genes for both categories
outlier_plot_combined <- ggplot(outlier_genes_combined, aes(x = Hugo_Symbol, y = Frequency, color = Category)) +
  geom_point(size = 3) +
  labs(title = "Outlier Mutation Frequencies", y = "Mutation Count", x = "Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print combined outlier plot
print(outlier_plot_combined)

# Save the combined outlier plot
ggsave("combined_outlier_mutations.png", outlier_plot_combined, width = 12, height = 8)

