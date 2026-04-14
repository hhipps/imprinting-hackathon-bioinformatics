# ---- 1. Load data ----
library(tidyverse)
library(readxl)

# the list of mammals
mammals <- read_excel("speciesTable.xlsx")
colnames(mammals)
head(mammals)

mammal_species <- mammals %>%
  filter(`NCBI class` == "Mammalia") %>% # includes marsupial and mammal in taxonomic group
  pull(`NCBI scientific name`) %>% # what is in the liver methyl matrix
  unique()
length(mammal_species)

# liver methylation matrix (genes x species)
methyl <- read_tsv("liver_methyl_matrix.txt") %>%
  rename(gene = `...1`) %>%
  filter(!is.na(gene)) %>%        # remove rows with no gene name
  mutate(gene = make.unique(gene)) %>%
  column_to_rownames(var = "gene")
head(rownames(methyl), 10)

# Read the imprinted genes file - irregular columns so read line by line
## could redo with info from geneimprint.com
lines <- readLines("ImpGenes Known.txt")

# Parse into a named list: species -> vector of gene symbols
imprinted_list <- lapply(lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  parts <- parts[parts != ""]  # remove empty fields
  list(species = parts[1], genes = parts[-1])
})

# Flatten to a two-column data frame: species | gene
imprinted_df <- do.call(rbind, lapply(imprinted_list, function(x) {
  if (length(x$genes) == 0) return(NULL)
  data.frame(species = x$species, gene = x$genes, stringsAsFactors = FALSE)
}))

# Get species-agnostic list of unique imprinted gene symbols
imprinted_genes_all <- unique(imprinted_df$gene)

# join with methyl matrix
# Add imprinting status based on row names (gene symbols)
methyl$is_imprinted <- rownames(methyl) %in% imprinted_genes_all
table(methyl$is_imprinted)
head(methyl$is_imprinted)


# filter methylation matrix by only mammals for initial checks
methyl_mammal <- methyl %>%
  rownames_to_column("gene") %>%
  select(gene, is_imprinted, any_of(mammal_species))

dim(methyl_mammal)
dim(methyl)
sum(colnames(methyl_mammal) %in% mammal_species) # how many species matched

# pivot mammal matrix to long format
methyl_mammal_long <- methyl_mammal %>%
  pivot_longer(-c(gene, is_imprinted), names_to = "species", values_to = "methylation") %>%
  filter(!is.na(methylation))
head(methyl_mammal_long) # compare to head(methyl_mammal)

# gene-level summary
gene_summary_mammal <- methyl_mammal_long %>%
  group_by(gene, is_imprinted) %>%
  summarise(mean_methyl = mean(methylation),
            median_methyl = median(methylation),
            var_methyl = var(methylation),
            n_species = n(),
            .groups = "drop")

#### Fig 1: DIST FROM 50% #####
# t-test on distance from 50%
gene_summary_mammal <- gene_summary_mammal %>%
  mutate(dist_from_50 = abs(mean_methyl - 50))

t.test(dist_from_50 ~ is_imprinted, data = gene_summary_mammal)

# histogram
ggplot(gene_summary_mammal, aes(x = mean_methyl, fill = is_imprinted)) +
  geom_histogram(binwidth = 5, position = "identity", alpha = 0.6) +
  scale_fill_manual(values = c("grey60", "steelblue"),
                    labels = c("Non-imprinted", "Imprinted")) +
  facet_wrap(~is_imprinted, labeller = label_both, scales = "free_y") +
  labs(x = "Mean Promoter Methylation in Mammals (%)",
       y = "Number of Genes", fill = NULL) +
  theme_minimal()


#### Fig. 2: VARIANCE ACROSS SPECIES #####
# t-test on variance across species
t.test(var_methyl ~ is_imprinted, data = gene_summary_mammal)

# violin plot to visualize variance difference
ggplot(gene_summary_mammal, aes(x = is_imprinted, y = var_methyl, fill = is_imprinted)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_manual(values = c("grey60", "steelblue"),
                    labels = c("Non-imprinted", "Imprinted")) +
  scale_x_discrete(labels = c("Non-imprinted", "Imprinted")) +
  labs(x = NULL, y = "Variance in Promoter Methylation Across Species",
       fill = NULL,
       title = "Methylation variance across mammalian species") +
  theme_minimal() +
  theme(legend.position = "none")

#### Fig. 3: LOGISTIC REGRESSION #####
install.packages("pROC")
library(pROC)

# 1. Prepare data
model_data <- gene_summary_mammal %>%
  filter(!is.na(var_methyl)) %>%  # remove genes with no variance (single species)
  mutate(is_imprinted = as.factor(is_imprinted))

# 2. Fit logistic regression
log_model <- glm(is_imprinted ~ mean_methyl + var_methyl,
                 data = model_data,
                 family = binomial)

summary(log_model)

# 3. Evaluate with ROC/AUC 
model_data$predicted_prob <- predict(log_model, type = "response")

roc_curve <- roc(model_data$is_imprinted, model_data$predicted_prob)
auc(roc_curve)

# plot ROC curve
ggroc(roc_curve, color = "steelblue", size = 1) +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed", color = "grey50") +
  annotate("text", x = 0.3, y = 0.2,
           label = paste0("AUC = ", round(auc(roc_curve), 3)),
           size = 5, color = "steelblue") +
  labs(title = "Logistic Regression: Predicting Imprinting Status",
       x = "Specificity", y = "Sensitivity") +
  theme_minimal()


