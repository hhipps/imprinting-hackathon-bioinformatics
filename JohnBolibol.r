library(tidyverse)
library(readxl)

animals <- read_excel("speciesTable.xlsx")

print(unique(animals$`NCBI class`))

mammals <- animals %>%
  filter(`NCBI class` == "Mammalia") %>%
  pull(`NCBI scientific name`) %>%
  unique()

head(mammals)
length(mammals)

methyl <- read_tsv("liver_methyl_matrix.txt") %>%
  rename(gene = `...1`) %>%
  filter(!(is.na(gene))) %>%
  mutate(gene = make.unique(gene)) %>%
  column_to_rownames(var = "gene")


# Imprinted gene list
lines <- readLines("ImpGenesKnown.txt")
imprinted_list <- lapply(lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  parts <- parts[parts != ""]
  list(species = parts[1], genes = parts[-1])
})
imprinted_df <- do.call(rbind, lapply(imprinted_list, function(x) {
  if (length(x$genes) == 0) return(NULL)
  data.frame(species = x$species, gene = x$genes, stringsAsFactors = FALSE)
}))
imprinted_genes_all <- unique(imprinted_df$gene)

methyl$is_imprinted <- rownames(methyl) %in% imprinted_genes_all

### Figure 1

methyl_mammal <- methyl %>%
  rownames_to_column("gene") %>%
  select(gene, is_imprinted, any_of(mammals))

methyl_mammal_long <- methyl_mammal %>%
  pivot_longer(-c(gene, is_imprinted), names_to = "species", values_to = "methylation") %>%
  filter(!(is.na(methylation)))

## summarize by gene
gene_summary_mammal <- methyl_mammal_long %>%
  group_by(gene, is_imprinted) %>%
  summarise(
    mean_methyl = mean(methylation),
    median_methyl = median(methylation),
    var_methyl = var(methylation),
    n_species = n(),
    .groups = "drop"
  ) %>%
  mutate(dist_from_50 = abs(mean_methyl - 50))

print(gene_summary_mammal)

fig1 = ggplot(gene_summary_mammal, aes(x = mean_methyl, fill = is_imprinted)) + 
  geom_histogram(binwidth = 5, position = "identity") +
  scale_fill_manual(values = c("black", "purple"), 
                    labels = c("Non-Imprinted", "Imprinted")) +
  facet_wrap(~is_imprinted, labeller = label_both, scales = "free_y") +
  labs(
    title = "Promoter Distribution in mammalian liver",
    subtitle = "Known imprinted vs non imprinted genes",
    x = "Mean Promoter Methylation Across Mammalian Species (%)",
    y = "Number of genes"
  )+
  theme_minimal()

  
fig1

head(methyl)
dim(methyl)