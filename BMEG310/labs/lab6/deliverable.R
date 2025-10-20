library(ISLR)

glm.fit <- glm(Direction ~ Lag1 + Lag2 + Lag3,
  data = Smarket,
  family = binomial
)

summary(glm.fit)

# %% Make predictions
glm.probs <- predict(glm.fit, type = "response")
glm.pred <- ifelse(glm.probs > 0.5, "Up", "Down")

confusion_matrix <- table(glm.pred, Smarket$Direction)
print("Confusion Matrix:")
print(confusion_matrix)

# Calculate classification rate (accuracy)
classification_rate <- mean(glm.pred == Smarket$Direction)
print(paste("Classification Rate:", round(classification_rate, 4)))
