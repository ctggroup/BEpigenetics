source("./enrichmentGO.R")
library(org.Hs.eg.db)

library(data.table)
egGO2ALLEGS<-getFromNamespace("org.Hs.egGO2ALLEGS","org.Hs.eg.db")
GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id","go_id","Ontology")]


load("./bmi_thinned/bmi_colnames_GS10k.RData")

C1_bmi <- fread(file = "./bmi_thinned/C1bmi.csv_thinned.csv")
C2_bmi <- fread(file = "./bmi_thinned/C2bmi.csv_thinned.csv")
C3_bmi <- fread(file = "./bmi_thinned/C3bmi.csv_thinned.csv")
C4_bmi <- fread(file = "./bmi_thinned/C4bmi.csv_thinned.csv")
print("enrichment BMI")
enrich_bmi<-enrichment(C1_bmi,C2_bmi,C3_bmi,C4_bmi,GeneID.PathID,probe_names=bmi_colnames,flatAnnon=flattenAnn("450k"),subsample=1:2) 
save(list="enrich_bmi",file="enrichmentGO_bmi.RData")


print("enrichment SMK")
C1_sm <- fread(file = "./sm_thinned/C1_sm.csv_thinned.csv")
C2_sm <- fread(file = "./sm_thinned/C2_sm.csv_thinned.csv")
C3_sm <- fread(file = "./sm_thinned/C3_sm.csv_thinned.csv")
C4_sm <- fread(file = "./sm_thinned/C4_sm.csv_thinned.csv")
enrich_smk<-enrichment(C1_sm,C2_sm,C3_sm,C4_sm,GeneID.PathID,probe_names=bmi_colnames,flatAnnon=flattenAnn("450k"),subsample=1:2) 
save(list="enrich_smk",file="enrichmentGO_smk.RData")
