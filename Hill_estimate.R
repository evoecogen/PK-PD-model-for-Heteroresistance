library(tidyverse)
library(drexplorer)
library(nlsMicrobio)

# Template for single antibiotic
# Dose-response with drug effects
# Input Dose response data
df_gr <- read_csv("Dose_response.csv")

# Convert it to data frame
df <- data.frame(x = df_gr$Concentration, y = df_gr$Growth)


# Set up Hill parameters and fitting Hill equation
m <- nls(y ~ (Gmax - ((Gmax-g)*(x/MIC)^H/((x/MIC)^H-g/Gmax))), data = df2, start = c(H = 1, g = -18, MIC = 1, Gmax = 1))
summary(m)

# Calculating R square
modelr::rsquare(m, df2)








