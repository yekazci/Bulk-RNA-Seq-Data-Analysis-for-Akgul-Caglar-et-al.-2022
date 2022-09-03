---
title: "Bulk RNA Seq data DEG Analysis"
author: "YEK"
date: "03 09 2022"
output: html_document
---


```{r setup, include = FALSE}
knitr::opts_chunk$set(eval = FALSE)
```

```{r}


 setwd("E:/TRANSCRIPTOME/1_Main_TRANSCRIPTOME/DRG_TRANSCRIPTOME_1_ANALYSIS/SALMON_DRG_NG_ANALYSIS_1")

 list.files()

 dir = getwd()

 salmon_drg_ng_samples = read.table(file.path(dir, "samples_drg_ng.txt"), header = TRUE)

 salmon_drg_ng_samples

 Samples

# 1 DRG_CS1_071217_S4_salmon_quant
# 2 DRG_CS2_071217_S5_salmon_quant
# 3 DRG_F_051017_S2_salmon_quant
# 4 DRG_F_121017_S1_salmon_quant
# 5 DRG_F_191017_S3_salmon_quant
# 6 NG_CS_071217_S8_salmon_quant
# 7 NG_CS_141217_S9_salmon_quant
# 8 NG_F_121017_S6_salmon_quant
# 9 G_F_191017_S7_salmon_quant

  sample_table_drg_ng = read.table("sample_table_drg_ng.txt", header = TRUE)

  sample_table_drg_ng

#   SampleNames     condition
# 1  DRG_CS1_S4      cont_drg
# 2  DRG_CS2_S5      cont_drg
# 3    DRG_F_S2    sample_drg
# 4    DRG_F_S1    sample_drg
# 5    DRG_F_S3    sample_drg
# 6    NG_CS_S8   cont_nodose
# 7    NG_CS_S9   cont_nodose
# 8     NG_F_S6 sample_nodose
# 9     NG_F_S7 sample_nodose

files <- file.path(dir,  "SALMON_DRG_RUN_1_QUANT_FILES", salmon_drg_ng_samples$Samples, "quant.sf")

   names(files) = paste0(sample_table_drg_ng$SampleNames)
 all(file.exists(files))

 # [1] TRUE

 files

# DRG_CS1_S4
# "E:/TRANSCRIPTOME/1_Main_TRANSCRIPTOME/DRG_TRANSCRIPTOME_1_ANALYSIS/SALMON_DRG_NG_ANALYSIS_1/SALMON_DRG_RUN_1_QUANT_FILES/DRG_CS1_071217_S4_salmon_quant/quant.sf"
# DRG_CS2_S5
# "E:/TRANSCRIPTOME/1_Main_TRANSCRIPTOME/DRG_TRANSCRIPTOME_1_ANALYSIS/SALMON_DRG_NG_ANALYSIS_1/SALMON_DRG_RUN_1_QUANT_FILES/DRG_CS2_071217_S5_salmon_quant/quant.sf"
# DRG_F_S2
# "E:/TRANSCRIPTOME/1_Main_TRANSCRIPTOME/DRG_TRANSCRIPTOME_1_ANALYSIS/SALMON_DRG_NG_ANALYSIS_1/SALMON_DRG_RUN_1_QUANT_FILES/DRG_F_051017_S2_salmon_quant/quant.sf"
# DRG_F_S1
# "E:/TRANSCRIPTOME/1_Main_TRANSCRIPTOME/DRG_TRANSCRIPTOME_1_ANALYSIS/SALMON_DRG_NG_ANALYSIS_1/SALMON_DRG_RUN_1_QUANT_FILES/DRG_F_121017_S1_salmon_quant/quant.sf"
# DRG_F_S3
# "E:/TRANSCRIPTOME/1_Main_TRANSCRIPTOME/DRG_TRANSCRIPTOME_1_ANALYSIS/SALMON_DRG_NG_ANALYSIS_1/SALMON_DRG_RUN_1_QUANT_FILES/DRG_F_191017_S3_salmon_quant/quant.sf"
# NG_CS_S8
# "E:/TRANSCRIPTOME/1_Main_TRANSCRIPTOME/DRG_TRANSCRIPTOME_1_ANALYSIS/SALMON_DRG_NG_ANALYSIS_1/SALMON_DRG_RUN_1_QUANT_FILES/NG_CS_071217_S8_salmon_quant/quant.sf"
# NG_CS_S9
# "E:/TRANSCRIPTOME/1_Main_TRANSCRIPTOME/DRG_TRANSCRIPTOME_1_ANALYSIS/SALMON_DRG_NG_ANALYSIS_1/SALMON_DRG_RUN_1_QUANT_FILES/NG_CS_141217_S9_salmon_quant/quant.sf"
# NG_F_S6
# "E:/TRANSCRIPTOME/1_Main_TRANSCRIPTOME/DRG_TRANSCRIPTOME_1_ANALYSIS/SALMON_DRG_NG_ANALYSIS_1/SALMON_DRG_RUN_1_QUANT_FILES/NG_F_121017_S6_salmon_quant/quant.sf"
# NG_F_S7
# "E:/TRANSCRIPTOME/1_Main_TRANSCRIPTOME/DRG_TRANSCRIPTOME_1_ANALYSIS/SALMON_DRG_NG_ANALYSIS_1/SALMON_DRG_RUN_1_QUANT_FILES/NG_F_191017_S7_salmon_quant/quant.sf"

 source("https://bioconductor.org/biocLite.R")

  biocLite("EnsDb.Mmusculus.v79")

#### different versions will be compared for tx2gene dataframes.

 ### make our own txdb for ensembl v91:

   library(ensembldb)

mouse_txdb_ens_v91 <- makeTxDbFromGFF("Mus_musculus.GRCm38.91.gtf/Mus_musculus.GRCm38.91.gtf",format="gtf")

# Import genomic features from the file as a GRanges object ... OK
# Prepare the 'metadata' data frame ... OK
# Make the TxDb object ... OK
# Warning message:
#   In .get_cds_IDX(type, phase) :
#   The "phase" metadata column contains non-NA values for features of type stop_codon. This information was ignored.

 keys_mouse_ens_v91 = keys(mouse_txdb_ens_v91, keytype = "GENEID")
 df_mouse_ens_91 <- select(mouse_txdb_ens_v91, keys = keys_mouse_ens_v91, keytype = "GENEID", columns = "TXNAME")

 # 'select()' returned 1:many mapping between keys and columns

 head(df_mouse_ens_91)

#  GENEID             TXNAME
# 1 ENSMUSG00000000001 ENSMUST00000000001
# 2 ENSMUSG00000000003 ENSMUST00000000003
# 3 ENSMUSG00000000003 ENSMUST00000114041
# 4 ENSMUSG00000000028 ENSMUST00000000028
# 5 ENSMUSG00000000028 ENSMUST00000096990
# 6 ENSMUSG00000000028 ENSMUST00000115585
 tx2gene_mouse_ens_91 = df_mouse_ens_91[, 2:1]  # tx ID, then gene ID
 head(tx2gene_mouse_ens_91)

#  TXNAME             GENEID
# 1 ENSMUST00000000001 ENSMUSG00000000001
# 2 ENSMUST00000000003 ENSMUSG00000000003
# 3 ENSMUST00000114041 ENSMUSG00000000003
# 4 ENSMUST00000000028 ENSMUSG00000000028
# 5 ENSMUST00000096990 ENSMUSG00000000028
# 6 ENSMUST00000115585 ENSMUSG00000000028

  length(tx2gene_mouse_ens_91$TXNAME)

  # [1] 133944

##### Other versions:

 library(TxDb.Mmusculus.UCSC.mm10.ensGene)

  keys_mouse_UCSC.mm10.ensGene = keys(TxDb.Mmusculus.UCSC.mm10.ensGene, keytype = "GENEID")

  df_mouse_UCSC.mm10.ensGene <- select(TxDb.Mmusculus.UCSC.mm10.ensGene, keys = keys_mouse_UCSC.mm10.ensGene, keytype = "GENEID", columns = "TXNAME")

  # 'select()' returned 1:many mapping between keys and columns

   head(df_mouse_UCSC.mm10.ensGene)

#    GENEID             TXNAME
# 1 ENSMUSG00000000001 ENSMUST00000000001
# 2 ENSMUSG00000000003 ENSMUST00000000003
# 3 ENSMUSG00000000003 ENSMUST00000114041
# 4 ENSMUSG00000000028 ENSMUST00000000028
# 5 ENSMUSG00000000028 ENSMUST00000096990
# 6 ENSMUSG00000000028 ENSMUST00000115585

 tx2gene_mouse_UCSC.mm10.ensGene = df_mouse_UCSC.mm10.ensGene[, 2:1]  # tx ID, then gene ID
 head(tx2gene_mouse_UCSC.mm10.ensGene)

#  TXNAME             GENEID
# 1 ENSMUST00000000001 ENSMUSG00000000001
# 2 ENSMUST00000000003 ENSMUSG00000000003
# 3 ENSMUST00000114041 ENSMUSG00000000003
# 4 ENSMUST00000000028 ENSMUSG00000000028
# 5 ENSMUST00000096990 ENSMUSG00000000028
# 6 ENSMUST00000115585 ENSMUSG00000000028

 length(tx2gene_mouse_UCSC.mm10.ensGene$TXNAME)

 # [1] 94647

  ### other version:

   library(EnsDb.Mmusculus.v79)
 txdb_mouse_ens_v79 = EnsDb.Mmusculus.v79
 keys_mouse_ens_v79 = keys(txdb_mouse_ens_v79, keytype = "GENEID")

   df_mouse_ens_v79 <- select(txdb_mouse_ens_v79, keys = keys_mouse_ens_v79, keytype = "GENEID", columns = "TXNAME")

   head(df_mouse_ens_v79)

#    GENEID             TXNAME
# 1 ENSMUSG00000000001 ENSMUST00000000001
# 2 ENSMUSG00000000003 ENSMUST00000000003
# 3 ENSMUSG00000000003 ENSMUST00000114041
# 4 ENSMUSG00000000028 ENSMUST00000000028
# 5 ENSMUSG00000000028 ENSMUST00000096990
# 6 ENSMUSG00000000028 ENSMUST00000115585

   tx2gene_mouse_ens_v79 = df_mouse_ens_v79[, 2:1]  # tx ID, then gene ID

   head(tx2gene_mouse_ens_v79)

#    TXNAME             GENEID
# 1 ENSMUST00000000001 ENSMUSG00000000001
# 2 ENSMUST00000000003 ENSMUSG00000000003
# 3 ENSMUST00000114041 ENSMUSG00000000003
# 4 ENSMUST00000000028 ENSMUSG00000000028
# 5 ENSMUST00000096990 ENSMUSG00000000028
# 6 ENSMUST00000115585 ENSMUSG00000000028

 length(tx2gene_mouse_ens_v79$TXNAME)

 # [1] 104129

 intersect_of_versions = intersect(tx2gene_mouse_ens_91$TXNAME, tx2gene_mouse_ens_v79$TXNAME)
 length(intersect_of_versions)

 # [1] 101977

 intersect_of_versions_2 = intersect(tx2gene_mouse_ens_91$TXNAME, tx2gene_mouse_UCSC.mm10.ensGene$TXNAME)
 length(intersect_of_versions_2)

 # [1] 91659

  #### Last version ens.v91 containes most transcript IDs in other previous versions.

   library(tximport)
 txi_mouse_ens_91_drg_nod = tximport(files, type = "salmon",ignoreTxVersion = TRUE ,tx2gene = tx2gene_mouse_ens_91)

#  reading in files with read_tsv
# 1 2 3 4 5 6 7 8 9
# transcripts missing from tx2gene: 74
# summarizing abundance
# summarizing counts
# summarizing length

   txi_mouse_ens_v79_drg_nod = tximport(files, type = "salmon",ignoreTxVersion = TRUE ,tx2gene = tx2gene_mouse_ens_v79)

#    reading in files with read_tsv
# 1 2 3 4 5 6 7 8 9
# transcripts missing from tx2gene: 24212
# summarizing abundance
# summarizing counts
# summarizing length

   txi_mouse_UCSC.mm10.ensGene_drg_nod = tximport(files, type = "salmon",ignoreTxVersion = TRUE ,tx2gene = tx2gene_mouse_UCSC.mm10.ensGene)

#    reading in files with read_tsv
# 1 2 3 4 5 6 7 8 9
# transcripts missing from tx2gene: 30599
# summarizing abundance
# summarizing counts
# summarizing length

   #### Number of missing transcripts are minimum in ensembl V.91 so we will go with it.

   sessionInfo()






```

