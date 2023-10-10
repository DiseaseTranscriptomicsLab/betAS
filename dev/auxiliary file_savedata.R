# Save data for r package

VT1_data_human <- readRDS("~/VersionControl/betAS/test/INCLUSION_LEVELS_FULL-hg19-98-v251.rds")
VT2_data_mouse <- readRDS("~/VersionControl/betAS/test/INCLUSION_LEVELS_FULL-mm10-8-v251.rds")
rMATS_data_mouse <- read.delim(file = "~/VersionControl/betAS/test/SE.MATS.JC.txt")
whippet_data_mouse <- readRDS("~/VersionControl/betAS/test/listdfs_WHippet.rds")

usethis::use_data(VT1_data_human)
usethis::use_data(VT2_data_mouse)
usethis::use_data(rMATS_data_mouse)
usethis::use_data(whippet_data_mouse)


VT1_metadata_human <- readRDS(file = "test/samplesTable.rds")
VT2_metadata_mouse <- readRDS(file = "test/samplesTable_Whippet.rds")
rMATS_metadata_mouse <- readRDS(file = "test/samplesTable_rMATS.rds")
whippet_metadata_mouse <- readRDS(file = "test/samplesTable_Whippet.rds")

usethis::use_data(VT1_metadata_human)
usethis::use_data(VT2_metadata_mouse)
usethis::use_data(rMATS_metadata_mouse)
usethis::use_data(whippet_metadata_mouse)


maxDevSimulationN100  <- readRDS("test/xintercepts_100incr_100cov_100trials.rds")
usethis::use_data(maxDevSimulationN100)
