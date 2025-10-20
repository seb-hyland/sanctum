data <- read.table("data.tsv", header = TRUE, sep = "\t")

print(dim(data))
print(colnames(data))

library(ggplot2)
ggplot(data, aes(x = Diagnosis, y = Age, fill = factor(Sex))) +
  geom_violin()
