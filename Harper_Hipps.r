install.packages("tidyverse")
install.packages("readxl")

library(tidyverse)
library(readxl)
library(pROC)

### get mammal list
animals <- read_excel("speciesTable.xlsx")

head(animals)

print(unique(animals$`NCBI class`))

mammals <- animals %>%
  filter(`NCBI class` == "Mammalia") %>%
  pull(`NCBI scientific name`) %>%
  unique()

head(mammals)
length(mammals)

### get methylation matrix
methyl <- read_tsv("liver_methyl_matrix.txt") %>%
  rename(gene = `...1`) %>%
  filter(!(is.na(gene))) %>%
  mutate(gene = make.unique(gene)) %>%
  column_to_rownames(var = "gene")
  
head(methyl)
dim(methyl)

### Imprinted gene list
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

### annotate with imprinting status
methyl$is_imprinted <- rownames(methyl) %in% imprinted_genes_all
View(methyl)

table(methyl$is_imprinted)
length(imprinted_genes_all)


### NEW Imprinted gene list
impNew <- read_csv("ImpGenesKnownNew.csv")
head(impNew)
length(unique(impNew$gene))
impSpecies <- unique(impNew$species)
head(methyl)

length(impSpecies)

matching_species <- intersect(impSpecies, colnames(methyl))
matching_species
length(matching_species)

missing_species <- setdiff(impSpecies, colnames(methyl))
missing_species

# remove any gene,species pairs that do not exist in the data
# Filter impNew to only rows where species is in your methyl dataset
impFiltered <- impNew %>% 
  filter(species %in% matching_species)

# Check how many genes remain
length(unique(impFiltered$gene))

# Compare to original
length(unique(impNew$gene))

############ FIGURE 1 #############

methyl_mammal <- methyl %>%
  rownames_to_column("gene") %>%
  select(gene, is_imprinted, any_of(mammals))

dim(methyl_mammal)
length(mammals)


methyl_mammal_long <- methyl_mammal %>%
  pivot_longer(-c(gene, is_imprinted), names_to = "species", values_to = "methylation") %>%
  filter(!(is.na(methylation)))

dim(methyl_mammal_long)

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

## make histogram
fig1 <- ggplot(gene_summary_mammal, aes(x = mean_methyl, fill = is_imprinted)) +
  geom_histogram(binwidth = 5, position = "identity") +
  scale_fill_manual(values = c("black", "purple"),
                    labels = c("Non-Imprinted", "Imprinted")) +
  facet_wrap(~is_imprinted, labeller = label_both, scales = "free_y") +
  labs(
    title = "Promoter Methylation Distribution in Mammalian Liver",
    subtitle = "Known imprinted vs non-imprinted genes",
    x = "Mean Promoter Methylation Across Mammalian Species (%)",
    y = "Number of genes"
  ) +
  theme_minimal()
fig1
ggsave("figures/fig1_methylation_distribution.png", fig1, width = 10, height = 5, dpi = 300)


## make histogram - dist from 50
fig05 <- ggplot(gene_summary_mammal, aes(x = dist_from_50, fill = is_imprinted)) +
  geom_histogram(binwidth = 5, position = "identity") +
  scale_fill_manual(values = c("black", "blue"),
                    labels = c("Non-Imprinted", "Imprinted")) +
  facet_wrap(~is_imprinted, labeller = label_both, scales = "free_y") +
  labs(
    title = "Promoter Methylation Distribution in Mammalian Liver",
    subtitle = "Known imprinted vs non-imprinted genes (dist from 50%)",
    x = "Distance of Mean Promoter Methylation from 50% (%)",
    y = "Number of genes"
  ) +
  theme_minimal()
fig05


## t test on mean methylation
cat("\n--- Figure 1: Mean Methylation t-test ---\n")
print(t.test(mean_methyl ~ is_imprinted, data = gene_summary_mammal))

## t test on diff from 50
cat("\n--- Figure 05: Distance from 50% t-test ---\n")
print(t.test(dist_from_50 ~ is_imprinted, data = gene_summary_mammal))

#is Macaca fascicularis in the mammals we have?

table(methyl_mammal_long$species == "Sus scrofa domesticus")

############ FIGURE 2 #############

fig2 <- ggplot(gene_summary_mammal, aes(x = is_imprinted, y = var_methyl, fill = is_imprinted)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  labs(
    title = "Methylation Variance Across Mammalian Species",
    x = NULL,
    y = "Variance in methylation"
  ) +
  scale_x_discrete(labels = c("Non-Imprinted", "Imprinted"))+
  theme_minimal() +
  theme(legend.position = "none")
fig2

ggsave("figures/fig2_methylation_variance_violin.png", fig2, width = 10, height = 5, dpi = 300)

# check how many genes are non-imprinted with variance > 2500
gene_summ_nonImp <- gene_summary_mammal %>%
  filter(is_imprinted == FALSE)
dim(gene_summ_nonImp)
sum(gene_summ_nonImp$var_methyl > 2500, na.rm = TRUE)

## t test on diff from 50
cat("\n--- Figure 2: Methylation variance t-test ---\n")
print(t.test(var_methyl ~ is_imprinted, data = gene_summary_mammal))

############ FIGURE 3 #############
## test cooks distance
varModel <- lm(gene_summary_mammal$var_methyl ~ gene_summary_mammal$is_imprinted)
cooks <- cooks.distance(varModel)

rowNum <- c(1:17370)
cookdf <- data.frame(rowNum, cooks)
cooks_plot <- ggplot(cookdf, aes(x=rowNum, y=cooks)) +
  geom_point() +
  labs(
    title = "Cook's Distance of Variance Linear Model",
    x = NULL,
    y = "Cook's Distance"
  ) +
  theme_minimal()
cooks_plot
ggsave("figures/fig3_cooks_dist.png", cooks_plot, width = 5, height = 5, dpi = 300)



############ FIGURE 4 #############
