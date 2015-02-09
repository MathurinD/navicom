#detach("package:RncMapping", unload=TRUE)
library(RncMapping)

conn = cBioConnect()

# Selection of the study id
studies = listStudies(conn, "leukemia")
st_id = studies$cancer_study_id[1]
ov_id = "ov_tcga_pub" # Published ovarian

visualizer = cBioNCviz(ov_id, genes_list="file://./genes_list", method="profiles")

saveInFiles(visualizer) # For NaviCell export
saveData(visualizer) # For easy sharing

