# %% PREP
df <- read.delim("ovarian.data", sep = ",", header = FALSE)
features <- c(
  "perimeter", "area", "smoothness", "symmetry", "concavity", paste(
    "protein", seq(1, 25),
    sep = ""
  )
)
names(df) <- c("cell_id", "diagnosis", features)
df$diagnosis <- as.factor(df$diagnosis)

features <- df[, 3:ncol(df)]

library(ggplot2)
library(nord)
library(ROCR)

head(features)



# %% Q1.1
pca <- prcomp(
  features,
  scale. = TRUE, center = TRUE
)
pca_summary <- summary(pca)
# Row 2 is variance
pc1_variance_prop <- pca_summary$importance[2, 1]
print(pc1_variance_prop)



# %% Q1.2
pca_cumvariance <- which(
  cumsum(pca_summary$importance[2, ]) >= 0.9
)[1]
print(pca_cumvariance)



# %% Q1.3
library(ggplot2)
coloring <- scale_color_manual(
  values = c("M" = "red", "B" = "blue"),
  labels = c("M" = "Malignant", "B" = "Benign"),
  name = "Diagnosis"
)

plot_data <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  diagnosis = df$diagnosis
)
ggplot(plot_data, aes(x = PC1, y = PC2, color = diagnosis)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "PC1/PC2 Features",
    x = "PC1",
    y = "PC2"
  ) +
  coloring



# %% Q1.4
feature_plot_data <- data.frame(
  area = df$area,
  concavity = df$concavity,
  diagnosis = df$diagnosis
)
ggplot(feature_plot_data, aes(x = area, y = concavity, color = diagnosis)) +
  geom_point(size = 3, alpha = 0.7) +
  labs(
    title = "Area vs Concavity Features",
    x = "Area",
    y = "Concavity"
  ) +
  coloring



# %% Q1.5
# TODO!



# %% Q1.6
pc_variance_prop <- pca_summary$importance[2, ]
pc_variance_df <- data.frame(
  PC = seq_along(pc_variance_prop),
  Value = pc_variance_prop
)
ggplot(pc_variance_df, aes(x = PC, y = Value, fill = as.factor(PC))) +
  geom_bar(stat = "identity") +
  labs(
    title = "Variance captured by each Principal Component",
    x = "Principal Component",
    y = "Percent variance",
  ) +
  scale_fill_nord("algoma_forest") +
  theme(legend.position = "none")




# %% Q2.1
metrics <- function(matrix) {
  tp <- matrix["M", "M"]
  tn <- matrix["B", "B"]
  fp <- matrix["M", "B"]
  fn <- matrix["B", "M"]

  accuracy <- (tp + tn) / sum(matrix)
  precision <- tp / (tp + fp)
  recall <- tp / (tp + fn)

  list(
    accuracy = accuracy,
    precision = precision,
    recall = recall
  )
}

print_metrics <- function(metrics) {
  for (name in names(metrics)) {
    val <- metrics[[name]]
    cat(name, ":", val, "\n")
  }
}

kmeans_analysis <- function(features, scaled = TRUE) {
  if (scaled) {
    features <- scale(features)
  }
  kmeans_result <- kmeans(features, centers = 2)
  clusters <- kmeans_result$cluster

  # First alignment
  pred_labels <- ifelse(clusters == 1, "M", "B")
  confusion_1 <- table(pred_labels, df$diagnosis)

  # Opposite alignment
  pred_labels_alt <- ifelse(clusters == 1, "B", "M")
  confusion_2 <- table(pred_labels_alt, df$diagnosis)

  # Higher diagonal sum = better alignment
  if (sum(diag(confusion_1)) >= sum(diag(confusion_2))) {
    confusion <- confusion_1
  } else {
    confusion <- confusion_2
  }

  # Metrics
  metrics <- metrics(confusion)

  list(
    confusion_1 = confusion_1,
    confusion_2 = confusion_2,
    metrics = metrics
  )
}

analysis <- kmeans_analysis(features)
kmeans_metrics <- analysis$metrics
print_metrics(kmeans_metrics)



# %% Q2.2
# Why identical???
replicates <- replicate(10, kmeans_analysis(features), simplify = FALSE)
replicate_metrics <- do.call(rbind, lapply(replicates, function(x) {
  as.data.frame(x$metrics)
}))
print(replicate_metrics)



# %% Q2.3
pca_features <- pca$x[, 1:5]
pca_analysis <- kmeans_analysis(pca_features, scaled = FALSE)
pca_metrics <- pca_analysis$metrics
print_metrics(pca_metrics)


# %% Q2.4
# TODO!




# %% Q3 PREP
training_set <- df[
  sample(nrow(df))
  [1:(nrow(df) / 2)],
]
testing_set <- df[
  sample(nrow(df))
  [(nrow(df) / 2):(nrow(df))],
]



# %% Q3.1
glm_fit <- glm(
  reformulate(
    names(training_set)[3:ncol(training_set)],
    response = "diagnosis"
  ),
  data = training_set,
  family = binomial
)


glm_probs <- predict(glm_fit, testing_set[3:ncol(testing_set)], type = "response")
glm_pred <- ifelse(glm_probs > 0.5, "M", "B")

confusion_matrix <- table(glm_pred, testing_set$diagnosis)
print("Confusion Matrix:")
print(confusion_matrix)

# Metrics
pred_metrics <- metrics(confusion_matrix)
print_metrics(pred_metrics)



# %% Q3.2
training_pca <- prcomp(
  training_set[, 3:ncol(training_set)],
  scale. = TRUE, center = TRUE
)
training_pca_df <- data.frame(
  diagnosis = training_set$diagnosis,
  training_pca$x[, 1:5]
)
testing_pca <- prcomp(
  testing_set[, 3:ncol(testing_set)],
  scale. = TRUE, center = TRUE
)
testing_pca_df <- data.frame(
  diagnosis = testing_set$diagnosis,
  testing_pca$x[, 1:5]
)

glm_fit_pca <- glm(
  reformulate(
    colnames(training_pca_df[, 2:ncol(training_pca_df)]),
    response = "diagnosis"
  ),
  data = training_pca_df,
  family = binomial
)


glm_probs_pca <- predict(glm_fit_pca, testing_pca_df[, 2:ncol(testing_pca_df)], type = "response")
glm_pred_pca <- ifelse(glm_probs_pca > 0.5, "M", "B")

confusion_matrix_pca <- table(glm_pred_pca, testing_set$diagnosis)
print("Confusion Matrix:")
print(confusion_matrix_pca)

# Metrics
pred_metrics <- metrics(confusion_matrix_pca)
print_metrics(pred_metrics)



# %% Q3.3
# TODO!



# %% Q3.4
# TODO!



# %% Q3.5
pred_prob <- predict(glm_fit, df[3:ncol(df)], type = "response")
predict <- prediction(pred_prob, df$diagnosis, label.ordering = c("B", "M"))
perform <- performance(predict, "tpr", "fpr")
plot(perform, colorize = TRUE)



# %% Q3.6
# %% Pt 1
lm_fit <- lm(
  reformulate(
    names(training_set)[3:ncol(training_set)],
    response = "diagnosis"
  ),
  data = training_set,
  family = binomial
)


lm_probs <- predict(lm_fit, testing_set[3:ncol(testing_set)], type = "response")
lm_pred <- ifelse(lm_probs > 0.5, "M", "B")

confusion_matrix <- table(lm_pred, testing_set$diagnosis)
print("Confusion Matrix:")
print(confusion_matrix)

# Metrics
pred_metrics <- metrics(confusion_matrix)
print_metrics(pred_metrics)


# %% Pt 2
lm_fit_pca <- lm(
  reformulate(
    colnames(training_pca_df[, 2:ncol(training_pca_df)]),
    response = "diagnosis"
  ),
  data = training_pca_df,
  family = binomial
)


lm_probs_pca <- predict(lm_fit_pca, testing_pca_df[, 2:ncol(testing_pca_df)], type = "response")
lm_pred_pca <- ifelse(lm_probs_pca > 0.5, "M", "B")

confusion_matrix_pca <- table(lm_pred_pca, testing_set$diagnosis)
print("Confusion Matrix:")
print(confusion_matrix_pca)

# Metrics
pred_metrics <- metrics(confusion_matrix_pca)
print_metrics(pred_metrics)
