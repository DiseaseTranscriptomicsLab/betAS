library(tidyverse)
library(highcharter)
library(hpackedbubble)

options(highcharter.theme = hc_theme_smpl(tooltip = list(valueDecimals = 2)))

testinggroupList <- list()
testinggroupList[[length(testinggroupList)+1]] <- list(name = "Control", samples = c("Control_1", "Control_2", "Control_3", "Control_4"), color = "#04A9FF")
testinggroupList[[length(testinggroupList)+1]] <- list(name = "Irrad", samples = c("Irrad_1", "Irrad_2", "Irrad_3", "Irrad_4"), color = "#DC143C")

samplelist    <- unlist(lapply(1:length(testinggroupList), function(x) testinggroupList[[x]]$samples))
groupList     <- unlist(lapply(1:length(testinggroupList), function(x) testinggroupList[[x]]$name))
sampleNumber  <- unlist(lapply(1:length(testinggroupList), function(x) length(testinggroupList[[x]]$samples)))


df <- data.frame("cat" = seq(1:length(samplelist)), "name" = rep(groupList, times = sampleNumber), "value" = samplelist)

# hc <- df %>%
#   hchart(
#     "pie", hcaes(x = name, y = y),
#     name = "Fruit consumption"
#   )


hpackedbubble(df$cat, df$name, df$value,
              title = "Groups",
              pointFormat = "<b>{point.name}:</b> {point.y}m CO<sub>2</sub>",
              dataLabelsFilter = 100,
              packedbubbleMinSize = "50%",
              packedbubbleMaxSize = "250%",
              theme = "sunset",
              packedbubbleZMin = 0,
              packedbubbleZmax = 10000, split = 0,
              gravitational = 0.02,
              parentNodeLimit = 1,
              dragBetweenSeries = 0,
              width = "100%")

data(gapminder, package = "gapminder")

gapminder <- gapminder %>%
  filter(year == max(year)) %>%
  select(country, pop, continent)

# hc <- hchart(gapminder, "packedbubble", hcaes(name = country, value = pop, group = continent))
# hc <- hchart(testinggroupList, "packedbubble", hcaes(name = name, value = samples, group = name))

q95 <- as.numeric(quantile(gapminder$pop, .95))

# hc %>%
#   hc_tooltip(
#     useHTML = TRUE,
#     pointFormat = "<b>{point.name}:</b> {point.value}"
#   ) %>%
#   hc_plotOptions(
#     packedbubble = list(
#       maxSize = "150%",
#       zMin = 0,
#       layoutAlgorithm = list(
#         gravitationalConstant =  0.05,
#         splitSeries =  TRUE, # TRUE to group points
#         seriesInteraction = TRUE,
#         dragBetweenSeries = TRUE,
#         parentNodeLimit = TRUE
#       ),
#       dataLabels = list(
#         enabled = TRUE,
#         format = "{point.name}",
#         filter = list(
#           property = "y",
#           operator = ">",
#           value = q95
#         ),
#         style = list(
#           color = "black",
#           textOutline = "none",
#           fontWeight = "normal"
#         )
#       )
#     )
#   )


#
# Highcharts.chart('container', {
#   chart: {
#     type: 'packedbubble',
#   },
#   series: [{
#     name: 'Coffee', // Coffee series
#     data: [{
#       // name property is used for the datalabel
#       // value property is used for the volume of the bubble
#       value: 12,
#       name: 'Bert'
#     }, {
#       value: 5,
#       name: 'John'
#     }, {
#       value: 10,
#       name: 'Sandra'
#     }, {
#       value: 7,
#       name: 'Cecile'
#     }]
#   }, {
#     name: 'Energy drinks', // Energy drinks series
#     data: [{
#       value: 10,
#       name: 'Tristan'
#     }]
#   }, {
#     name: 'Tea', // Tea series
#     data: [5, 6, 8, {
#       value: 10,
#       name: 'Mustapha',
#       color: 'pink'
#     }]
#   }]
# });
