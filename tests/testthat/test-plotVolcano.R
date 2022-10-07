test_that("multiplication works", {
  testTable <- betAS:::testTable
  maxDevSimulationN100  <- betAS:::maxDevSimulationN100
  psiTable   <- testTable$PSI
  qualTable  <- testTable$Qual
  testGroups <- list()
  testGroups[[length(testGroups)+1]] <- list(name = "GroupA", samples = c("ERR2598266", "ERR2598267", "ERR2598268"), color = "#FF9AA2")
  testGroups[[length(testGroups)+1]] <- list(name = "GroupB", samples = c("ERR2598270", "ERR2598307", "ERR2598351"), color = "#FFB7B2")
  groupA    <- "GroupA"
  groupB    <- "GroupB"
  samplesA   <- testGroups[[1]]$samples
  samplesB   <- testGroups[[2]]$samples
  colsGroupA    <- convertCols(psiTable, samplesA)
  colsGroupB    <- convertCols(psiTable, samplesB)
  volcanoTable <- prepareTableVolcano(psitable = psiTable, qualtable = qualTable, npoints = 500, colsA = colsGroupA, colsB = colsGroupB, labA = groupA, labB = groupB, basalColor = "#89C0AE", interestColor = "#E69A9C", maxDevTable = maxDevSimulationN100)

  expect_s3_class(volcanoTable, "data.frame")
  expect_equal(volcanoTable[ , 1:ncol(psiTable)], cbind(psiTable))
  expect_equal(colnames(volcanoTable), c(colnames(psiTable), "Pdiff", "betasPsiA", "betasPsiB", "deltapsi"))
  expect_false(any(is.na(volcanoTable)))
  expect_equal(volcanoTable[["deltapsi"]], volcanoTable[["betasPsiB"]] - volcanoTable[["betasPsiA"]])

  plot <- plotVolcano(betasTable = volcanoTable, labA = groupA, labB = groupB, basalColor = "#89C0AE", interestColor = "#E69A9C")
  expect_s3_class(plot, "ggplot")
  expect_equal(plot$data, volcanoTable)
  expect_equal(length(plot$layers), 3)

})
