#change working directory to where "Heloderma_transcript_parse.sh" was run,
#recommended running in R Studio, but can be run as "Rscript Heloderma_dosage_violins.R" in an environment such as,
#"mamba create -n vioplot R r-vioplot"
#"conda activate vioplot"

#set working directory, e.g.
#setwd("~/Heloderma_violins")

#Final violin plot uncorrected data
#install.packages("vioplot")
library("vioplot")


#Figure 3 plot
svg("rplot1.svg") #<-- Rscript adaptation

par(mfrow=c(2,1))
par(mar=c(3,4,2.0,2.0))

gila_auto <- read.delim("gila2.refbased.stringtie_compiled_transcripts_auto.txt", header = T)
gallus_auto <- read.delim("z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_auto.txt", header = T)
gila_sex <- read.delim("gila2_galgal5_stringtie_transcripts.txt", header = T)
gallus_28 <- read.delim("galgal5_gila2_stringtie_transcripts.txt", header = T)
gallus_sex  <- read.delim("z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_sex.txt", header = T)

F_M_list <- list(Gal_auto = log2(gallus_auto$f_m_ratio_refbased), Gal_28 = log2(gallus_28$f_m_ratio_refbased), 
                 Gal_Z = log2(gallus_sex$f_m_ratio_refbased), Helo_auto = log2(gila_auto$f_m_ratio_refbased), 
                 Helo_Z = log2(gila_sex$f_m_ratio_refbased))
F_M_means <- as.data.frame(sapply(F_M_list,mean))
F_M_medians <- as.data.frame(sapply(F_M_list,median))

vioplot(F_M_list, ylim = range(-5:5), 
        ylab = "log2(F/M ratio of ref-based transcript expression)", 
        main = "Gallus + Heloderma F/M expression", names = F)
points(F_M_means, pch = 19, col = "grey", cex = 1)
legend("topright", pch = c(21, 19), col = c("black", "grey"),
       bg = "white", legend = c("Median", "Mean"), cex = 1)


wilcox.test(F_M_list$Helo_auto, F_M_list$Helo_Z, paired = F, alternative = "two.sided")
#p-value = p-value = 1.5e-13
wilcox.test(F_M_list$Gal_auto, F_M_list$Gal_Z, paired = F, alternative = "two.sided")
#p-value < 2.2e-16
wilcox.test(F_M_list$Gal_auto, F_M_list$Gal_28, paired = F, alternative = "two.sided")
#p-value = 0.9301


male_female_list <- list(Helo_auto_M = log2(gila_auto$male_mean_refbased), Helo_auto_F = log2(gila_auto$female_mean_refbased), 
                         Helo_Z_M = log2(gila_sex$male_mean_refbased), Helo_Z_F = log2(gila_sex$female_mean_refbased), 
                         Gal_auto_M = log2(gallus_auto$male_mean_refbased), Gal_auto_F = log2(gallus_auto$female_mean_refbased), 
                         Gal_28_M = log2(gallus_28$male_mean_refbased),  Gal_28_F = log2(gallus_28$female_mean_refbased), 
                         Gal_Z_M = log2(gallus_sex$male_mean_refbased),  Gal_Z_F = log2(gallus_sex$female_mean_refbased))

vioplot(male_female_list$Gal_auto_M, male_female_list$Gal_28_M, male_female_list$Gal_Z_M, male_female_list$Helo_auto_M, male_female_list$Helo_Z_M, side = "left", plotCentre = "line", col = "grey40", 
        ylab = "log2(ref-based transcript expression)", main = "Gallus + Heloderma, F + M expression", names = F, ylim = range(-15:15))
vioplot(male_female_list$Gal_auto_F, male_female_list$Gal_28_F, male_female_list$Gal_Z_F, male_female_list$Helo_auto_F, male_female_list$Helo_Z_F, side = "right", plotCentre = "line", col = "grey80", add = T)
legend("topright", col = c("grey40", "grey80"), pch = 15,
       bg = "white", legend = c("Male", "Female"), cex = 1.5)
mtext("Gallus auto                                                              Gallus 28                                                          Gallus Z                                                            Heloderma auto                                                            Heloderma Z", side=1, line=1)


male_female_sum <- as.data.frame(sapply(male_female_list,mean))
male_female_sum$median <- cbind(sapply(male_female_list,median))


wilcox.test(male_female_list$Gal_Z_M, male_female_list$Gal_Z_F, paired = F, alternative = "two.sided")
#p-value = 0.006835

wilcox.test(male_female_list$Gal_auto_M, male_female_list$Gal_28_M, paired = F, alternative = "two.sided")
#p-value = 0.01592

wilcox.test(male_female_list$Gal_auto_F, male_female_list$Gal_28_F, paired = F, alternative = "two.sided")
#p-value = 0.01819

wilcox.test(male_female_list$Helo_Z_M, male_female_list$Helo_Z_F, paired = F, alternative = "two.sided")
#p-value = 0.1077

wilcox.test(male_female_list$Helo_auto_M, male_female_list$Helo_auto_F, paired = F, alternative = "two.sided")
#p-value = 0.4268

wilcox.test(male_female_list$Helo_Z_M, male_female_list$Helo_auto_M, paired = F, alternative = "two.sided")
#p-value = 0.09872

wilcox.test(male_female_list$Helo_auto_F, male_female_list$Helo_Z_F, paired = F, alternative = "two.sided")
#p-value = 0.7751

dev.off() #<-- Rscript adaptation



#Supplemental Figure 1 plot
svg("rplot2.svg")  #<-- Rscript adaptation


#box plot with Heloderma genes clustered by Gallus chromosomes
#setwd("E:/Drive/Current_Manuscripts/Heloderma/final")

par(mfrow=c(2,1))
par(mar=c(3,4,2.0,2.0))

chr1  <- read.delim("gila2.refbased.stringtie_chr1.txt", header = T)
chr2  <- read.delim("gila2.refbased.stringtie_chr2.txt", header = T)
chr3  <- read.delim("gila2.refbased.stringtie_chr3.txt", header = T)
chr4  <- read.delim("gila2.refbased.stringtie_chr4.txt", header = T)
chr5  <- read.delim("gila2.refbased.stringtie_chr5.txt", header = T)
chr6  <- read.delim("gila2.refbased.stringtie_chr6.txt", header = T)
chr7  <- read.delim("gila2.refbased.stringtie_chr7.txt", header = T)
chr8  <- read.delim("gila2.refbased.stringtie_chr8.txt", header = T)
chr9  <- read.delim("gila2.refbased.stringtie_chr9.txt", header = T)
chr10 <- read.delim("gila2.refbased.stringtie_chr10.txt", header = T)
chr11 <- read.delim("gila2.refbased.stringtie_chr11.txt", header = T)
chr12 <- read.delim("gila2.refbased.stringtie_chr12.txt", header = T)
chr13 <- read.delim("gila2.refbased.stringtie_chr13.txt", header = T)
chr14 <- read.delim("gila2.refbased.stringtie_chr14.txt", header = T)
chr15 <- read.delim("gila2.refbased.stringtie_chr15.txt", header = T)
#chr16 <- read.delim("gila2.refbased.stringtie_chr16.txt", header = T)
#chr17 <- read.delim("gila2.refbased.stringtie_chr17.txt", header = T)
chr18 <- read.delim("gila2.refbased.stringtie_chr18.txt", header = T)
chr19 <- read.delim("gila2.refbased.stringtie_chr19.txt", header = T)
chr20 <- read.delim("gila2.refbased.stringtie_chr20.txt", header = T)
chr21 <- read.delim("gila2.refbased.stringtie_chr21.txt", header = T)
chr22 <- read.delim("gila2.refbased.stringtie_chr22.txt", header = T)
chr23 <- read.delim("gila2.refbased.stringtie_chr23.txt", header = T)
chr24 <- read.delim("gila2.refbased.stringtie_chr24.txt", header = T)
chr25 <- read.delim("gila2.refbased.stringtie_chr25.txt", header = T)
chr26 <- read.delim("gila2.refbased.stringtie_chr26.txt", header = T)
chr27 <- read.delim("gila2.refbased.stringtie_chr27.txt", header = T)
chr28 <- read.delim("gila2.refbased.stringtie_chr28.txt", header = T)
chr33 <- read.delim("gila2.refbased.stringtie_chr33.txt", header = T)
chrZ  <- read.delim("gila2.refbased.stringtie_chrZ.txt", header = T)

chr_list <- list(	
  chr1  = log2(chr1$f_m_ratio_refbased),
  chr2  = log2(chr2$f_m_ratio_refbased),
  chr3  = log2(chr3$f_m_ratio_refbased),
  chr4  = log2(chr4$f_m_ratio_refbased),
  chr5  = log2(chr5$f_m_ratio_refbased),
  chr6  = log2(chr6$f_m_ratio_refbased),
  chr7  = log2(chr7$f_m_ratio_refbased),
  chr8  = log2(chr8$f_m_ratio_refbased),
  chr9  = log2(chr9$f_m_ratio_refbased),
  chr10 = log2(chr10$f_m_ratio_refbased),
  chr11 = log2(chr11$f_m_ratio_refbased),
  chr12 = log2(chr12$f_m_ratio_refbased),
  chr13 = log2(chr13$f_m_ratio_refbased),
  chr14 = log2(chr14$f_m_ratio_refbased),
  chr15 = log2(chr15$f_m_ratio_refbased),
#  chr16 = log2(chr16$f_m_ratio_refbased),
#  chr17 = log2(chr17$f_m_ratio_refbased),
  chr18 = log2(chr18$f_m_ratio_refbased),
  chr19 = log2(chr19$f_m_ratio_refbased),
  chr20 = log2(chr20$f_m_ratio_refbased),
  chr21 = log2(chr21$f_m_ratio_refbased),
  chr22 = log2(chr22$f_m_ratio_refbased),
  chr23 = log2(chr23$f_m_ratio_refbased),
  chr24 = log2(chr24$f_m_ratio_refbased),
  chr25 = log2(chr25$f_m_ratio_refbased),
  chr26 = log2(chr26$f_m_ratio_refbased),
  chr27 = log2(chr27$f_m_ratio_refbased),
  chr28 = log2(chr28$f_m_ratio_refbased),
  chr33 = log2(chr33$f_m_ratio_refbased),
  chrZ  = log2(chrZ$f_m_ratio_refbased))

library("vioplot")

vioplot(chr_list, ylim = range(-5:5), 
        ylab = "log2(F/M FPKM ratio)", main = "Heloderma F:M expression", names = F)
mtext("chr1   chr2   chr3   chr4   chr5   chr6   chr7   chr8   chr9   chr10   chr11   chr12   chr13   chr14   chr15   chr18   chr19   chr20   chr21   chr22   chr23   chr24   chr25   chr26   chr27   chr28   chr33   chrZ", side=1, line=1)



Gallus_transcripts <- read.delim("z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_all.txt", header = T)

chr_list <- list(	
  chr1  = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr1", 6]),
  chr2  = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr2", 6]),
  chr3  = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr3", 6]),
  chr4  = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr4", 6]),
  chr5  = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr5", 6]),
  chr6  = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr6", 6]),
  chr7  = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr7", 6]),
  chr8  = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr8", 6]),
  chr9  = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr9", 6]),
  chr10 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr10", 6]),
  chr11 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr11", 6]),
  chr12 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr12", 6]),
  chr13 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr13", 6]),
  chr14 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr14", 6]),
  chr15 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr15", 6]),
# chr16 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr16", 6]),
# chr17 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr16", 6]),
  chr18 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr18", 6]),
  chr19 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr19", 6]),
  chr20 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr20", 6]),
  chr21 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr21", 6]),
  chr22 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr22", 6]),
  chr23 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr23", 6]),
  chr24 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr24", 6]),
  chr25 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr25", 6]),
  chr26 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr26", 6]),
  chr27 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr27", 6]),
  chr28 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr28", 6]),
  chr33 = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chr33", 6]),
  chrZ  = log2(Gallus_transcripts[Gallus_transcripts$scaffold == "chrZ", 6]))

vioplot(chr_list, ylim = range(-5:5), 
        ylab = "log2(F/M FPKM ratio)", main = "Gallus F:M expression", names = F)
mtext("chr1   chr2   chr3   chr4   chr5   chr6   chr7   chr8   chr9   chr10   chr11   chr12   chr13   chr14   chr15   chr18   chr19   chr20   chr21   chr22   chr23   chr24   chr25   chr26   chr27   chr28   chr33   chrZ", side=1, line=1)

dev.off()  #<-- Rscript adaptation
