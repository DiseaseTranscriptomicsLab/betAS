# Testing VT functions ----------------------------------------------------


filterVT <- list()

testTable <- read.table(gzfile("/home/mariana/betAS/test/INCLUSION_LEVELS_FULL-Hsa32-hg19_to_test.tab.gz"), sep="\t", header=TRUE, quote="")
# maxDevSimulationN100 <- readRDS(url("http://imm.medicina.ulisboa.pt/group/distrans/SharedFiles/Mariana/Splicing&SenescenceFLEX/xintercepts_100incr_100cov_100trials.R"))
maxDevSimulationN100  <- readRDS("test/xintercepts_100incr_100cov_100trials.rds")

#mimic filtering from function
incTable <- testTable
qual_cols <- grep("-Q", colnames(incTable))
psi_cols  <- qual_cols-1

psiVAST  <- incTable[,c(1:6,psi_cols)]
qualVAST <- incTable[,c(1:6,qual_cols)]

nrow(psiVAST)
originalColN <- ncol(psiVAST)

# *** to be made REACTIVE!!! ***
# Consider exon skipping events only
psiVAST <- psiVAST[which(psiVAST$COMPLEX == "S"),]
nrow(psiVAST)

# Remove events containing at least one NA
psiVAST$AnyNA    <- apply(psiVAST, 1, anyNA)
psiVAST          <- psiVAST[which(psiVAST$AnyNA == FALSE),]
psiVAST          <- psiVAST[,-c(ncol(psiVAST))]
qualVAST         <- qualVAST[match(psiVAST$EVENT, qualVAST$EVENT),]
nrow(psiVAST)

# Filter out events from "ANN" module out
# annot             <- which(psiVAST$COMPLEX == "ANN")
# psiVAST           <- psiVAST[-c(annot),]
# qualVAST          <- qualVAST[match(psiVAST$EVENT, qualVAST$EVENT),]
# nrow(psiVAST)

# Calculate coverage/balance (use .Q columns)
qualVAST$AllminVLOW  <- apply(qualVAST[,grep(".Q", colnames(qualVAST))], 1, VT_all_minVLOW_tags)
qualVAST             <- qualVAST[which(qualVAST$AllminVLOW == TRUE),]
psiVAST              <- psiVAST[match(qualVAST$EVENT, psiVAST$EVENT),]

# Remove columns added
psiVAST               <- psiVAST[,c(1:originalColN)]
qualVAST              <- qualVAST[,c(1:originalColN)]

filterVT[[1]] <- psiVAST
filterVT[[2]] <- qualVAST
filterVT[[3]] <- table(psiVAST$COMPLEX)
filterVT[[4]] <- colnames(psiVAST)[-c(1:6)]

names(filterVT) <- c("PSI", "Qual", "EventsPerType", "Samples")
return(filterVT)


# Simulate filtering
testTable <- read.table(gzfile("/home/mariana/betAS/test/INCLUSION_LEVELS_FULL-Hsa32-hg19_to_test.tab.gz"), sep="\t", header=TRUE, quote="")
colnames(testTable) <- gsub(x = colnames(testTable), pattern = "Sample_IMR90_", replacement = "")
samples <- colnames(testTable)[grep(x = colnames(testTable), pattern = "-Q")-1]
cat(paste0("Initial events: ", nrow(testTable)))

test <- filterVastTools(testTable)
cat(paste0("Filtered events: ", nrow(test$PSI)))

altPsi <- alternativeVastTools(test, minPsi = 5, maxPsi = 95)
cat(paste0("Alternative events: ", nrow(altPsi$PSI)))

# Prepare "Big picture plot"
#
# bigPicturePlot <- function(table){
#
#   transf <- reshape2::melt(table[,-c(1,3:6)], id.vars = c("EVENT"))
#
#   plot <- ggplot(data = transf,
#          aes(x = value,
#              group = variable,
#              color = variable,
#              fill = variable)) +
#     geom_density(show.legend = FALSE,
#                  alpha = 0.5) +
#     scale_x_continuous(breaks = seq(0,100, 25), limits = c(0,100)) +
#     facet_wrap(. ~ variable, ncol = 4) +
#     theme_minimal()
#
#   return(plot)
# }

transf <- reshape2::melt(altPsi$PSI[,-c(1,3:6)], id.vars = c("EVENT"))
colors <- colorRampPalette(c("#55CC9D","#F2969A"))(32)

ggplot(data = transf,
       aes(x = value,
           group = variable,
           color = variable,
           fill = variable)) +
  geom_density(show.legend = FALSE,
               alpha = 0.8,
               size = 2) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = seq(0,100, 25), limits = c(0,100)) +
  facet_wrap(. ~ variable, ncol = 4) +
  theme_betAS()

table(transf$variable)
samplesGroupsList <- list(`Control` = c("Control_1", "Control_2", "Control_3", "Control_4"),
                          `Doxorubicin` = c("Doxo_1", "Doxo_2", "Doxo_3", "Doxo_4"),
                          `Etoposide` = c("Eto_1", "Eto_2", "Eto_3", "Eto_4"),
                          `Palbociclib` = c("Palbo_1", "Palbo_2", "Palbo_3", "Palbo_4"),
                          `Irradiation` = c("Irrad_1", "Irrad_2", "Irrad_3", "Irrad_4"),
                          `Ras+` = c("iRas_pos_rel_1", "iRas_pos_rel_2", "iRas_pos_rel_3", "iRas_pos_rel_4"),
                          `Ras-` = c("iRas_neg_rel_1", "iRas_neg_rel_2", "iRas_neg_rel_3", "iRas_neg_rel_4"),
                          `Quiescence` = c("Quiesc_1", "Quiesc_2", "Quiesc_3", "Quiesc_4"))
exeventnames <- c("HsaEX0007927", "HsaEX0032264", "HsaEX0039848", "HsaEX0029465", "HsaEX0026102", "HsaEX0056290", "HsaEX0035084", "HsaEX0065983", "HsaEX0036532", "HsaEX0049206")


#Prepare individual betAS plot
maxDevSimulationN100 <- readRDS(url("http://imm.medicina.ulisboa.pt/group/distrans/SharedFiles/Mariana/Splicing&SenescenceFLEX/xintercepts_100incr_100cov_100trials.R"))

#inputs for global densities (betAS)
eventID <- "HsaEX0055568"
npoints <- 500
psitable <- test$PSI
qualtable <- test$Qual
colsA <- c(7:10)
colsB <- c(27:30)
labA <- "Control"
labB <- "Irrad"
colorA <- c("#04A9FF")
colorB <- c("#B3B3B3")

plotIndividualDensities(eventID = "HsaEX0055568", npoints = 500, psitable = altPsi$PSI, qualtable = altPsi$Qual, colsA = c(7:10), colsB = c(27:30), labA = "Control", labB = "Irrad", colorA = c("#B3B3B3"), colorB =  c("#04A9FF"))

#inputs for global densities (betAS)
eventID <- "HsaEX0055568"
npoints <- 500
psitable <- test$PSI
qualtable <- test$Qual
testinggroupList <- list()
testinggroupList[[length(testinggroupList)+1]] <- list(name = "Control", samples = c("Control_1", "Control_2", "Control_3", "Control_4"), color = "#04A9FF")
testinggroupList[[length(testinggroupList)+1]] <- list(name = "Irrad", samples = c("Irrad_1", "Irrad_2", "Irrad_3", "Irrad_4"), color = "#DC143C")

plotIndividualDensitiesList(eventID = "HsaEX0055568", npoints = 500, psitable = altPsi$PSI, qualtable = altPsi$Qual, groupList = testinggroupList)


# Run betAS for all events in table
psitable
qualtable
npoints
colsA
colsB
themeColors <- c("#89C0AE", "#E69A9C", "#76CAA0", "#EE805B", "#F7CF77", "#81C1D3")

plotVolcano(psitable = test$PSI, qualtable = test$Qual, npoints = 500, colsA = c(7:10), colsB = c(27:30), basalColor = themeColors[1], interestColor = themeColors[4], maxDevTable = maxDevSimulationN100)

# Testing group betAS
testList <- list()
testList[[length(testList)+1]] <- list(name = "Control", samples = c("Control_1", "Control_2", "Control_3", "Control_4"), color = "#04A9FF")
testList[[length(testList)+1]] <- list(name = "Test1", samples = c("Test_1", "Test_2"), color = "#DC143C")

description <- character()

for(i in 1:length(testList)){

  concatenate <- function(elem){paste0(c("<ul>", "<h1>", elem$name, "</h1>", "<li>", as.character(elem$samples), "</li>", "</ul>"))}

  cat(unlist(sapply(testList, concatenate)))

  description <- paste0(description, "<ul>", "<h1>", testList[[i]]$name, "</h1>", "<li>", as.character(testList[[i]]$samples), "</li>", "</ul>")

}

description <- paste0("<h1>", description, "</h1>")

concatenate <- function(elem){paste0(c("<ul>", "<h1>", elem$name, "</h1>", "<li>", paste(elem$samples, sep = " PLEASE "), "</li>", "</ul>"))}

# Convert groups list into table
unlist(testList)

table <- c()
for(i in 1:length(testList)){

  group <- c("Name" = testList[[i]]$name, "Samples" = paste0(testList[[i]]$samples, collapse = " "), "Color" = testList[[i]]$color)
  table <- rbind(table, group)

}





# Implement F-stat in the volcano plot
psitable = test$PSI
qualtable = test$Qual
npoints = 500
colsA = c(7:10)
colsB = c(27:30)
labA = "Control"
labB = "Irrad"
basalColor = themeColors[1]
interestColor = themeColors[4]
maxDevTable = maxDevSimulationN100
groupList = testinggroupList
# plotVolcanoFstat <- function(psitable, qualtable, npoints, colsA, colsB, labA, labB, basalColor, interestColor, maxDevTable){

colsA    <- convertCols(psitable, colsA)
samplesA <- names(colsA)

colsB    <- convertCols(psitable, colsB)
samplesB <- names(colsB)

#Group betAS (A)
groupAbetAS <- lapply(1:nrow(qualtable),
                      function(x)
                        individualBetas_nofitting_incr(table = qualtable[x,],
                                                       cols = colsA,
                                                       indpoints = npoints,
                                                       maxdevRefTable = maxDevTable))

#Group betAS (B)
groupBbetAS <- lapply(1:nrow(qualtable),
                      function(x)
                        individualBetas_nofitting_incr(table = qualtable[x,],
                                                       cols = colsB,
                                                       indpoints = npoints,
                                                       maxdevRefTable = maxDevTable))

#Assign names based on event IDs
names(groupAbetAS) <- qualtable$EVENT
names(groupBbetAS) <- qualtable$EVENT


# WORKING ON FDR APPROACH
testEvent <- "HsaEX0055568"

#Group betAS (A)
groupAbetAS <- individualBetas_nofitting_incr(table = qualtable[grep(testEvent, qualtable$EVENT),],
                                                       cols = colsA,
                                                       indpoints = npoints,
                                                       maxdevRefTable = maxDevTable)
#Group betAS (B)
groupBbetAS <- individualBetas_nofitting_incr(table = qualtable[grep(testEvent, qualtable$EVENT),],
                                                       cols = colsB,
                                                       indpoints = npoints,
                                                       maxdevRefTable = maxDevTable)

indBetasA <- groupAbetAS
indBetasB <- groupBbetAS
groupsAB  <- c(labA, labB)


estimateFDR <- function(indBetasA, indBetasB, groupsAB){

  nsim <- 1000
  jntBetas <- list()

  # Get predefined objects A:
  artifA <- indBetasA$BetaPoints
  artiflabelsA <- names(indBetasA$BetaPoints)
  psisA <- indBetasA$PSI
  nA <- length(indBetasA$PSI)
  labelA <- groupsAB[1]
  covA <- indBetasA$inc + indBetasA$exc
  medianA <- indBetasA$MedianBeta

  # Get predefined objects B:
  artifB <- indBetasB$BetaPoints
  artiflabelsB <- names(indBetasB$BetaPoints)
  psisB <- indBetasB$PSI
  nB <- length(indBetasB$PSI)
  labelB <- groupsAB[2]
  covB <- indBetasB$inc + indBetasB$exc
  medianB <- indBetasB$MedianBeta

  # Under the null hypothesis of no difference in PSI across groups (all samples come from the same distribution)
  originalMedian <- median(c(indBetasA$PSI, indBetasB$PSI))
  foundPsi <- medianB - medianA

  simDeltaPsi <- c()
  for(i in 1:nsim){

    simulatedA <- c()
    for(sample in unique(artiflabelsA)){

      simReads <- simulate_reads(cov = covA[which(names(covA) == sample)], psi = originalMedian)
      simulatedA <- c(simulatedA, rbeta(500, shape1 = simReads$inc, shape2 = simReads$exc))

    }

    simulatedB <- c()
    for(sample in unique(artiflabelsB)){

      simReads <- simulate_reads(cov = covB[which(names(covB) == sample)], psi = originalMedian)
      simulatedB <- c(simulatedB, rbeta(500, shape1 = simReads$inc, shape2 = simReads$exc))

    }

    deltaPsi <- median(simulatedB) - median(simulatedA)
    simDeltaPsi <- c(simDeltaPsi, deltaPsi)

  }

  fdr <- length(abs(simDeltaPsi)<= abs(foundPsi))/nsim

  plot <- ggplot() +
    geom_density(aes(x = simDeltaPsi)) +
    geom_vline(xintercept = medianB - medianA,
               color = "#DC143C",
               size = 2) +
    xlab(expression(Delta*PSI[simulated])) +
    # annotate(geom = "text") +
    ylab("") +
    theme_betAS()

  jntBetas[[1]] <- simDeltaPsi
  jntBetas[[2]] <- foundPsi
  jntBetas[[3]] <- fdr
  jntBetas[[4]] <- plot

  names(jntBetas) <- c("Simulated dPSI", "Found dPSI", "FDR", "Plot")
  return(jntBetas)

}

astring <- colnames(psitable)[7]
bstring <- colnames(psitable)[8]
longest_string <- function(s){return(s[which.max(nchar(s))])}
adist(x = a, y = b)
agrep(pattern = a, x = b)

nchar(a)
nchar(b)
## get all forward substrings of 'b'
sb <- stri_sub(b, 1, 1:nchar(b))
## extract them from 'a' if they exist
sstr <- na.omit(stri_extract_all_coll(a, sb, simplify=TRUE))
## match the longest one
sstr[which.max(nchar(sstr))]





agrep("lasy", "1 lazy 2")
agrep("lasy", c(" 1 lazy 2", "1 lasy 2"), max.distance = list(sub = 0))
agrep("laysy", c("1 lazy", "1", "1 LAZY"), max.distance = 2)
agrep("laysy", c("1 lazy", "1", "1 LAZY"), max.distance = 2, value = TRUE)
agrep("laysy", c("1 lazy", "1", "1 LAZY"), max.distance = 2, ignore.case = TRUE)

exnames <- colnames(psitable)[-c(1:6)]

table = psitable
findGroupsVast <- function(table){



}




# Testing Whippet functions ----------------------------------------------------

folder_Whippet <- "/mnt/scratch/home/mariana/Projects/betAS/Whippet/"
files <- paste0("/mnt/scratch/home/mariana/Projects/betAS/Whippet/",list.files(folder_Whippet, pattern = "\\.psi\\.gz$"))
listfiles <- lapply(files,fread)
names(listfiles) <- sapply(files, function(file) gsub("\\..*","",gsub(".*/","",file)))

#saveRDS(object = listfiles, file = "test/listdfs_WHippet.rds")
whip1 <- getWhippet(listfiles)

filtwhip1 <- filterWhippet(whip1, names(whip1$EventsPerType))

altfiltwhip1 <- alternativeWhippet(filtwhip1,1,99)


metadatawhippet <- readRDS("test/samplesTable_rMATS.rds")
metadatawhippet <- metadatawhippet[,-1]
metadatawhippet
#saveRDS(object = metadatawhippet, file = "test/samplesTable_Whippet.rds")


#testonesamp <- readRDS("test/INCLUSION_LEVELS_FULL-hg19-98-v251.rds")
#testonesamp <- testonesamp[,c(1:8)]
#saveRDS(testonesamp, "test/testVT_onesamp.tab")


# Testing rMATS functions ----------------------------------------------------
testTable <- read.delim(file = "test/SE.MATS.JC.txt")
testTable <- getrMATS(testTable)
testTable <- filterrMATS(testTable)
testTable <- alternativerMATS(testTable, minPsi=1, maxPsi=100)

maxDevSimulationN100  <- readRDS("test/xintercepts_100incr_100cov_100trials.rds")
row=1
columns=c(7,8,9,10)
npoints=500

individualBetas_nofitting_incr(table = testTable$Qual[row,],
                               cols = columns,
                               indpoints = npoints,
                               maxdevRefTable = maxDevSimulationN100)



testTable

RM_all_minReads()

quals <- testTable$Qual[1,7:ncol(testTable$Qual)]

# # Saving metadata for HS dataset (rMATS)
# aux <- as.data.frame(fread("~/p16_ARF/Hernandez_Segura/metadata.txt"))
# aux <- aux[match(unique(aux$`Comment[ENA_RUN]`), aux$`Comment[ENA_RUN]`),]
# aux$sampleID <- paste0("G1_S",1:nrow(aux))
# aux <- aux[,c("sampleID","Comment[ENA_RUN]","Factor Value[cell type]","Factor Value[irradiate]")]
#
# colnames(aux) <- c("sampleID","Run","CellType","Irradiate")
# row.names(aux) <- 1:nrow(aux)
# #three samples per condition
# subsetsamplesidx <- c(1,2,3,31,39,50,55,63,74,4,9,18,5,6,7,32,33,34,56,57,58,8,10,11,38,40,41,62,64,65,12,13,14,45,46,47,69,70,71)
# subsetsamplesidx <- subsetsamplesidx[order(subsetsamplesidx)]
# aux <- aux[subsetsamplesidx,]
# table(aux$CellType,aux$Irradiate)
#
# saveRDS(object = aux, file = "test/samplesTable_rMATS_HS.rds")
#
#
#
#
# folder_CardosoMoreira <- "/mnt/nmorais-nfs/Backup_Unit/rsilva/RitaS_backup/Data/HernandezSegura/HernandezSeguraFASTQ/"
#
# files <- paste0("/mnt/nmorais-nfs/Backup_Unit/rsilva/RitaS_backup/Data/HernandezSegura/HernandezSeguraFASTQ/",list.files(folder_CardosoMoreira, pattern = "\\.fastq\\.gz$"))
# files_1 <- files[grep("_1.fastq.gz",files)]
# files_2 <- files[grep("_2.fastq.gz",files)]
#
# files_CardosoMoreira <- paste0(files_1,":",files_2)
#
#
# files_CardosoMoreira <- unique(grep(paste(aux$Run,collapse="|"),
#                         files_CardosoMoreira, value=TRUE))
#
#
# files_CardosoMoreira <- paste(files_CardosoMoreira, sep="", collapse=",")
#
# write.table(files_CardosoMoreira, file = "/mnt/scratch/home/rsilva/betAS_tests/rMATS/HernandezSegura/group1.txt", sep = "\t",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)


