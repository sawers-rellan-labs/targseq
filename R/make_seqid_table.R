library(dplyr)

#data from
# https://docs.google.com/spreadsheets/d/1gsZ017XvS_xLZkXzBmKT7e4DL5aiEld8AcxfeJ2KNwc/edit?gid=1345415053#gid=1345415053&fvid=1451101345
line_info <- read_csv("~/Desktop/FINAL_SELECTION.csv")
sample_annotation <- read_csv("~/Desktop/sample_annotation.csv")

seqid_label <-  line_info %>%
  select(sample, donor_accession = accession, taxa,introgression=hpc1, label_1 =label) %>%
  mutate(maizegdb_prefix = gsub("\\..*","",donor_accession, perl = TRUE),
         seqid =paste("hpc1",sample, sep="_"), 
         seqid_sample=sample,
         sample_n = gsub("S","",sample, fixed = TRUE),
         donor_sufix = gsub(".*\\.","",donor_accession, perl = TRUE),
         ancestry_preffix  = c("R","D")[introgression+1],
         ancestry = c("Recurrrent","Donor")[introgression+1]) %>%
  select(seqid, seqid_sample, sample_n, maizegdb_prefix,taxa, donor_accession, taxa, donor_sufix,ancestry,ancestry_preffix,label_1)



seqid_label <- rbind( seqid_label %>%
  mutate(
    label_2 =paste(ancestry_preffix,maizegdb_prefix, taxa,donor_sufix, seqid_sample, sep = "_"),
    label_3 =paste(maizegdb_prefix, taxa,donor_sufix, seqid_sample, sep = "_")
    #   paste( ancestry,maizegdb_prefix, taxa,donor_sufix,sample,sep="_",colapse=TRUE),
  ) %>% select(seqid,founder_ancestry = ancestry,ancestry_preffix, starts_with("label")),
  # add information for maizegdb genomes, they are not samples
data.frame(seqid	="hpc1_B73", founder_ancestry="Recurrent",ancestry_preffix="R",
           label_1 ="B73", label_2="R_Zm_B73",	label_3= "Zm_B73"),
data.frame(seqid	="hpc1_TIL18", founder_ancestry="Donor",ancestry_preffix="D",
           label_1 ="TIL18", label_2="D_Zx_TIL18",	label_3= "Zx_TIL18")
)

write.csv(seqid_label,"~/Desktop/seqid_label_1.csv",quote=FALSE,row.names = FALSE)
