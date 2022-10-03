install.packages("gprofiler2")
install.packages('jsonlite', version = '1.8.1') #DO NOT COMPILE JSONLITE

library(gprofiler2)

gostres <- gost(query = c("RMRP", "DBET", "H2BC18", "H4C5", "LGR4-AS1","FSTL1","BCYRN1","H4C8","LSM8","PNPLA3","ANKRD34A","COTL1","HEATR3","GREM1","RPPH1",
                          "ANKRD30B","TUBA1B","TMSB10","ACTG1","MT-ND3","HSPA8","MT-CYB","MT-ND4L","FASTK","TNPO2","CTTN-DT","GPM6B","GAPDH","VEZF1",
                          "MT-CO2","MT-ND4","TAOK1","MT-CO3","PKM","H2AC8","H2BC9","HNRNPH1","ANKMY1","RUFY2","SPAG9","LNPK","GGA3","PLEKHG2","RSPH1","RNF133",
                          "MTF1","ST6GALNAC2","WTIP","KLHL25"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = FALSE, exclude_iea = TRUE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = c("GO:CC", "GO:BP", "GO:MF"), as_short_link = FALSE)

p <- gostplot(gostres, capped = FALSE, interactive = FALSE)
p
pp <- publish_gostplot(p, highlight_terms = c("GO:0098803", "GO:0005746", "GO:0070469", "GO:0098800", "GO:0043505", "GO:0061638", "GO:0034506", "GO:0000786", "GO:0098798", "GO:0070069", "GO:1990204", "GO:0005747", "GO:0030964","GO:0045271"), 
                       width = NA, height = NA, filename = NULL )

