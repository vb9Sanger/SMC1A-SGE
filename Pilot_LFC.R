if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

exon3.1 <- readRDS("sgRNA_HG_16_1190757649_library_independent_deseq2_results.rds")
View(exon3.1[["Library Sequence LFC Simplified Table"]])
SMCA1_exon3_sg1 <- exon3.1[["Library Sequence LFC Simplified Table"]]
write.csv(SMCA1_exon3_sg1, file='/Users/vb9/Downloads/SMC1A_exon3_sg1.csv')
View(exon3.1[["Library Sequence LFC Main Table"]])

exon3.2 <- readRDS("sgRNA_HG_16_1190757657_library_independent_deseq2_results.rds")
SMC1A_exon3_sg2 <- exon3.2[["Library Sequence LFC Simplified Table"]]
write.csv(SMC1A_exon3_sg2, file='/Users/vb9/Downloads/SMC1A_exon3_sg2.csv')

exon4.1 <- readRDS("sgRNA_HG_168_1190757455_library_independent_deseq2_results.rds")
SMC1A_exon4_sg1 <- exon4.1[["Library Sequence LFC Simplified Table"]]
write.csv(SMC1A_exon4_sg1, file='/Users/vb9/Downloads/SMC1A_exon4_sg1.csv')

exon4.2 <- readRDS("sgRNA_HG_168_1190757473_library_independent_deseq2_results.rds")
SMC1A_exon4_sg2 <- exon4.2[["Library Sequence LFC Simplified Table"]]
write.csv(SMC1A_exon4_sg2, file='/Users/vb9/Downloads/SMC1A_exon4_sg2.csv')

