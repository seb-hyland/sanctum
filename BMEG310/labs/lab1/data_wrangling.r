# %%%% Introduction to data wrangling and functions in R
# %%
age <- c(15, 22, 45, 52, 73, 81)
age[5]
age[-5]
age[c(3, 5, 6)]
age[1:4]

# %%
age > 50
age > 50 | age < 18
age[age > 50 | age < 18]
idx <- age > 50 | age < 18
age[idx]

# %% Dataframes
v <- read.csv("mouse_exp_design.csv")
View(v)
genotype <- v["genotype"]
genotype[sample6]
