library(readxl)
data <- read_excel("dataset.xls", sheet = "Hoja3")
model <- lm(height ~ age + playtime, data = data)

plot(model$fitted.values, model$residuals,
  xlab = "Fitted Values",
  ylab = "Residuals",
  main = "Residual Plot: Height vs Predicted Height"
)
