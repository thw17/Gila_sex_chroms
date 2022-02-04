#change working directory to where "Heloderma_transcript_parse.sh" was run,
#recommended running in R Studio, but can be run as "Rscript Heloderma_dosage_violins.R" in an environment such as,
#"mamba create -n vioplot R r-vioplot"
#"conda activate vioplot"

#set working directory, e.g.
#setwd("~/Heloderma_violins")

#Final violin plot uncorrected data
#install.packages("vioplot")
library("vioplot")

#svg("rplot.svg") #<-- Rscript adaptation

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
legend("topleft", pch = c(21, 19), col = c("black", "grey"),
       bg = "white", legend = c("Median", "Mean"), cex = 1)

abline(h = median(F_M_list$Helo_auto))
abline(h = median(F_M_list$Gal_auto), col = "grey")
names <- c("Heloderma F/M Autosomal Median", "Gallus F/M Autosomal Median")
colors <- c("black", "grey")
legend("topright", legend=names,
       col=colors, lty=1, cex=0.9,
       text.font=4)

wilcox.test(F_M_list$Helo_auto, F_M_list$Helo_Z, paired = F, alternative = "two.sided")
#p-value = 0.00000000000015
wilcox.test(F_M_list$Gal_auto, F_M_list$Gal_Z, paired = F, alternative = "two.sided")
#p-value < 0.00000000000000022



male_female_list <- list(Helo_auto_M = log2(gila_auto$male_mean_refbased), Helo_auto_F = log2(gila_auto$female_mean_refbased), 
                         Helo_Z_M = log2(gila_sex$male_mean_refbased), Helo_Z_F = log2(gila_sex$female_mean_refbased), 
                         Gal_auto_M = log2(gallus_auto$male_mean_refbased), Gal_auto_F = log2(gallus_auto$female_mean_refbased), 
                         Gal_28_M = log2(gallus_28$male_mean_refbased),  Gal_28_F = log2(gallus_28$female_mean_refbased), 
                         Gal_Z_M = log2(gallus_sex$male_mean_refbased),  Gal_Z_F = log2(gallus_sex$female_mean_refbased))

vioplot(male_female_list$Gal_auto_M, male_female_list$Gal_28_M, male_female_list$Gal_Z_M, male_female_list$Helo_auto_M, male_female_list$Helo_Z_M, side = "left", plotCentre = "line", col = "grey40", 
        ylab = "log2(ref-based transcript expression)", main = "Gallus + Heloderma, F + M expression", names = F, ylim = range(-15:15))
vioplot(male_female_list$Gal_auto_F, male_female_list$Gal_28_F, male_female_list$Gal_Z_F, male_female_list$Helo_auto_F, male_female_list$Helo_Z_F, side = "right", plotCentre = "line", col = "grey80", add = T)
legend("topleft", col = c("grey40", "grey80"), pch = 15,
       bg = "white", legend = c("Male", "Female"), cex = 1.5)
mtext("Gallus auto                                                              Gallus 28                                                          Gallus Z                                                            Heloderma auto                                                            Heloderma Z", side=1, line=1)
t.test(male_female_list$Gal_Z_M, male_female_list$Gal_Z_F)
#p-value = 0.01734
t.test(male_female_list$Helo_Z_M, male_female_list$Helo_Z_F)
#p-value = 0.5309
abline(h = 0)
abline(h = mean(male_female_list$Helo_Z_M), col = "grey80")
abline(h = mean(male_female_list$Helo_Z_F), col = "grey40")
names <- c("Heloderma M Z Mean", "Heloderma F Z Mean", "'0' gridline")
colors <- c("grey40", "grey80", "black")
legend("topright", legend=names,
       col=colors, lty=1, cex=0.9,
       text.font=4)

male_female_sum <- as.data.frame(sapply(male_female_list,mean))
male_female_sum$median <- cbind(sapply(male_female_list,median))


wilcox.test(male_female_list$Gal_Z_M, male_female_list$Gal_Z_F, paired = F, alternative = "two.sided")
#p-value = 0.006835
wilcox.test(male_female_list$Helo_Z_M, male_female_list$Helo_Z_F, paired = F, alternative = "two.sided")
#p-value = 0.1077
wilcox.test(male_female_list$Helo_Z_M, male_female_list$Helo_auto_M, paired = F, alternative = "two.sided")
#p-value = 0.09872
wilcox.test(male_female_list$Helo_auto_F, male_female_list$Helo_Z_F, paired = F, alternative = "two.sided")
#p-value = 0.7751

#dev.off() #<-- Rscript adaptation































#quick analysis with "corrected" data:
#library("vioplot")
#par(mfrow=c(2,1))
#par(mar=c(3,4,2.0,2.0))
#
##gallus_auto_uncorr <- read.delim("z_ortho_filtered-gila2_galgal5-galgal5_refbased.stringtie_compiled_per_transcript_separate_individuals.txt")
##gallus_28_uncorr   <- read.delim("z_ortho_filtered-gila2_galgal5-galgal5_refbased.stringtie_compiled_per_transcript_separate_individuals.txt")
##gallus_sex_uncorr  <- read.delim("z_ortho_filtered-gila2_galgal5-galgal5_refbased.stringtie_compiled_per_transcript_separate_individuals.txt")
#
#gallus_auto_corr <- read.delim("z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_auto.txt")
#gallus_28_corr <- read.delim("galgal5_gila2_stringtie_transcripts_corrected.txt")
#gallus_sex_corr <- read.delim("z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_sex.txt")
#gila_auto_corr <- read.delim("corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_auto.txt")
#gila_sex_corr  <- read.delim("gila2_galgal5_stringtie_transcripts_corrected.txt")
#
#Gal_corr_list <- list(Gal_auto = log2(gallus_auto_corr$corrected_f_m_ratio), Gal_28 = log2(gallus_28_corr$corrected_f_m_ratio), 
#                 Gal_Z = log2(gallus_sex_corr$corrected_f_m_ratio), Helo_auto = log2(gila_auto_corr$corrected_f_m_ratio), 
#                 Helo_Z = log2(gila_sex_corr$corrected_f_m_ratio))
#Gal_corr_means <- as.data.frame(sapply(Gal_corr_list,mean))
#Gal_corr_medians <- as.data.frame(sapply(Gal_corr_list,median))
#
#vioplot(Gal_corr_list, ylim = range(-5:5), 
#        ylab = "log2(F/M ratio of ref-based transcript expression)", 
#        main = "Corrected Gallus + Heloderma F/M expression", names = F)
#points(Gal_corr_means, pch = 19, col = "grey", cex = 1)
#legend("topleft", pch = c(21, 19), col = c("black", "grey"),
#       bg = "white", legend = c("Median", "Mean"), cex = 1)
#abline(h = median(Gal_corr_list$Helo_auto))
#abline(h = median(Gal_corr_list$Gal_auto), col = "grey")
#names <- c("Heloderma F/M Autosomal Median", "Gallus F/M Autosomal Median")
#colors <- c("black", "grey")
#legend("topright", legend=names,
#       col=colors, lty=1, cex=0.9,
#       text.font=4)
#
#
#male_female_corr_list <- list(Helo_auto_M = log2(gila_auto_corr$corrected_male_mean), Helo_auto_F = log2(gila_auto_corr$corrected_female_mean), 
#                         Helo_Z_M = log2(gila_sex_corr$corrected_male_mean), Helo_Z_F = log2(gila_sex_corr$corrected_female_mean), 
#                         Gal_auto_M = log2(gallus_auto_corr$corrected_male_mean), Gal_auto_F = log2(gallus_auto_corr$corrected_female_mean), 
#                         Gal_28_M = log2(gallus_28_corr$corrected_male_mean),  Gal_28_F = log2(gallus_28_corr$corrected_female_mean), 
#                         Gal_Z_M = log2(gallus_sex_corr$corrected_male_mean),  Gal_Z_F = log2(gallus_sex_corr$corrected_female_mean))
#male_female_corr_sum <- as.data.frame(sapply(male_female_corr_list,mean))
#male_female_corr_sum$median <- cbind(sapply(male_female_corr_list,median))
#
#vioplot(male_female_corr_list$Gal_auto_M, male_female_corr_list$Gal_28_M, male_female_corr_list$Gal_Z_M, male_female_corr_list$Helo_auto_M, male_female_corr_list$Helo_Z_M, side = "left", plotCentre = "line", col = "grey40", 
#        ylab = "log2(ref-based transcript expression)", main = "Corrected Gallus + Heloderma, F + M expression", names = F, ylim = range(-15:15))
#vioplot(male_female_corr_list$Gal_auto_F, male_female_corr_list$Gal_28_F, male_female_corr_list$Gal_Z_F, male_female_corr_list$Helo_auto_F, male_female_corr_list$Helo_Z_F, side = "right", plotCentre = "line", col = "grey80", add = T)
#legend("topleft", col = c("grey40", "grey80"), pch = 15,
#       bg = "white", legend = c("Male", "Female"), cex = 1.5)
#mtext("Gallus auto                                                              Gallus 28                                                          Gallus Z                                                            Heloderma auto                                                            Heloderma Z", side=1, line=1)
#abline(h = 0)
#abline(h = mean(male_female_corr_list$Helo_Z_M), col = "grey80")
#abline(h = mean(male_female_corr_list$Helo_Z_F), col = "grey40")
#names <- c("Heloderma M Z Mean", "Heloderma F Z Mean", "'0' gridline")
#colors <- c("grey40", "grey80", "black")
#legend("topright", legend=names,
#       col=colors, lty=1, cex=0.9,
#       text.font=4)
#