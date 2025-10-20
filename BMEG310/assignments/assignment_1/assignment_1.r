data <- read.csv(file = "rna_seq.csv", header = TRUE)


# %% Q1
median_expression <- median(data[["expression"]], na.rm = TRUE)
above_median_expression_subset <- subset(data, expression > median_expression)

unique_genes <- length(unique(above_median_expression_subset[["gene"]]))
print(unique_genes)


# %% Q2
female_infected_subset <- subset(
  data,
  infection == "InfluenzaA" & sex == "Female"
)
unique_mice <- length(unique(female_infected_subset[["mouse"]]))
print(unique_mice)


# %% Q3
expressed_subset <- subset(data, expression > 0)
male_unique_subset <- unique(
  subset(expressed_subset, sex == "Male")[["gene"]]
)
female_unique_subset <- unique(
  subset(expressed_subset, sex == "Female")[["gene"]]
)
common_genes <- length(
  intersect(male_unique_subset, female_unique_subset)
)
print(common_genes)


# %% Q4
chromosome_one_subset <- subset(data, chromosome_name == 1 & expression > 0)
list_per_mouse <- split(chromosome_one_subset, chromosome_one_subset[["mouse"]])
lowest_expression_genes <- lapply(list_per_mouse, function(mouse_df) {
  mouse_df[which.min(mouse_df[["expression"]]), c("mouse", "gene")]
})
recombined_df <- do.call(rbind, lowest_expression_genes)
colnames(recombined_df) <- c("mouse", "lowest_gene")
print(recombined_df, row.names = FALSE)


# %% Q5
get_top_3_chromosomes <- function(input_df) {
  list_per_chromosome <- split(input_df, input_df[["chromosome_name"]])
  list_per_chromosome_top3 <- lapply(
    list_per_chromosome, function(chromosome_df) {
      chromosome_df[
        order(
          chromosome_df[["expression"]],
          decreasing = TRUE
        )[1:3], c("chromosome_name", "gene", "expression")
      ]
    }
  )
  recombined_chromosome_df <- do.call(rbind, list_per_chromosome_top3)
  recombined_chromosome_df
}
result <- get_top_3_chromosomes(data)
print(result, row.names = FALSE)


# %% Q6
library(ggplot2)
library(nord)
sex_factor <- as.factor(data$sex)
infection_factor <- as.factor(data$infection)
ggplot(
  data,
  aes(x = sex_factor, y = expression, fill = infection_factor, theme = "dark")
) +
  theme_bw() +
  labs(
    title = "Gene expression levels based on sex and infection",
    y = "Expression (" ~ log[10] ~ "scale)",
    x = "Sex",
    fill = "Infection"
  ) +
  geom_violin(scale = "count") +
  scale_y_continuous(transform = "log10") +
  scale_fill_nord("algoma_forest")
