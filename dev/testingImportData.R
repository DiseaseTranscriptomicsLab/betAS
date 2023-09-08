# Testing import data

# Examples provided by betAS ---------------------------------------------------
# rMATS
test_rMATS <- getDataset(pathTables=NULL, tool="rMATS") # input should be a path to a *MATS.JC.txt
test_rMATS_PSI <- getEvents(test_rMATS, tool="rMATS")
test_rMATS_filter <- filterEvents(test_rMATS_PSI, types=NULL, N=10 )
test_rMATS_alt <- alternativeEvents(test_rMATS_filter, minPsi=0, maxPsi=100 )

# vast-tools
test_VT <- getDataset(pathTables=NULL, tool="vast-tools") # input should be a path to a *INCLUSION_LEVELS_FULL*.tab file
test_VT_PSI <- getEvents(test_VT, tool="vast-tools")
test_VT_filter <- filterEvents(test_VT_PSI, types=NULL, N=10 )
test_VT_alt <- alternativeEvents(test_VT_filter, minPsi=0, maxPsi=100 )

# whippet
test_whippet <- getDataset(pathTables=NULL, tool="whippet") # input should be a list of paths to *.psi.gz files
test_whippet_PSI <- getEvents(test_whippet, tool="whippet")
test_whippet_filter <- filterEvents(test_whippet_PSI, types=NULL, N=10 )
test_whippet_alt <- alternativeEvents(test_whippet_filter, minPsi=0, maxPsi=100 )

# User-input data --------------------------------------------------------------
# rMATS
test_rMATS_UI <- getDataset(pathTables="test/SE.MATS.JC.txt", tool="rMATS") # input should be a path to a *MATS.JC.txt
test_rMATS_PSI_UI <- getEvents(test_rMATS_UI, tool="rMATS")
test_rMATS_filter_UI <- filterEvents(test_rMATS_PSI_UI, types=NULL, N=10 )
test_rMATS_alt_UI <- alternativeEvents(test_rMATS_filter_UI, minPsi=0, maxPsi=100  )
#
# # vast-tools
# test_VT <- getDataset(pathTables=NULL, tool="vast-tools") # input should be a path to a *INCLUSION_LEVELS_FULL*.tab file
# test_VT_PSI <- getEvents(test_VT, tool="vast-tools")
# test_VT_filter <- filterEvents(test_VT_PSI, types=NULL, N=10)
# test_VT_alt <- alternativeEvents(test_VT_filter, minPsi=0, maxPsi=100  )
#
# # whippet
# test_whippet <- getDataset(pathTables=NULL, tool="whippet") # input should be a list of paths to *.psi.gz files
# test_whippet_PSI <- getEvents(test_whippet, tool="whippet")
# test_whippet_filter <- filterEvents(test_whippet_PSI, types=NULL, N=10 )
# test_whippet_alt <- alternativeEvents(test_whippet_filter, minPsi=0, maxPsi=100 )


