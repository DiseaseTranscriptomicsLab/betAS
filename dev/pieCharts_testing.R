# Load required R packages
library(dplyr)
library(highcharter)

# from https://www.datanovia.com/en/lessons/highchart-interactive-pie-chart-and-alternatives-in-r/
# Set highcharter options
options(highcharter.theme = hc_theme_smpl_tailored(tooltip = list(valueDecimals = 2)))

df <- data.frame(
  x = c(0, 1, 2, 3, 4),
  y = c(10, 19.4, 21.1, 14.4, 6.4),
  name = as.factor(c("grape", "olive", "guava", "nut", "pear"))
)
df

hc <- df %>%
  hchart(
    "pie", hcaes(x = name, y = y),
    name = "Fruit consumption"
  )

hc

table(psitable$COMPLEX)

preparePieForVastToolsCOMPLEX <- function(table){

  # Set highcharter options
  options(highcharter.theme = hc_theme_smpl_tailored(tooltip = list(valueDecimals = 2)))

  df <- as.data.frame(table(table$COMPLEX))
  df <- df[,c(2,1)]
  df <- cbind(seq(0, nrow(df)-1), df)
  colnames(df) <- c("x", "y", "name")

  hc <- df %>%
    hchart(
      "pie", hcaes(x = name, y = y),
      name = "Fruit consumption"
    )

  return(hc)

}
