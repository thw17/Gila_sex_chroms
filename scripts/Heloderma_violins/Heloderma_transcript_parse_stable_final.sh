#!bin/bash
#3 files needed in the current directory for this code to run:
#"z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt.gz"
#"gila2_galgal5_gff_comparison.txt.gz"
#"gila2.refbased.stringtie_compiled_per_transcript.txt.gz"
#"chr.txt.gz"
#"gila_chicken_gff_comparison_FULL.txt.gz"

gunzip *.gz

#make header for parsed files
head -1 z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt > z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_chr.txt
head -1 z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt > z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_auto.txt
head -1 z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt > z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_sex.txt
#add filtered data to header files, filtering out non-chromosomes, the W, the mtDNA, etc and transcripts that aren't expressed in either sex
grep "chr" z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt | grep -v "_" | grep -v "[[:space:]]0.0[[:space:]]" | grep -v "M" | grep -v "W" | grep -v "LGE64" >> z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_chr.txt
grep "chr" z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt | grep -v "_" | grep -v "[[:space:]]0.0[[:space:]]" | grep -v "28" | grep -v "M" | grep -v "W" | grep -v "Z" | grep -v "LGE64" >> z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_auto.txt
grep "chrZ" z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt | grep -v "_" | grep -v "[[:space:]]0.0[[:space:]]" >> z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_sex.txt
#remove "chr" before the chromosome names
sed -i 's/chr//g' z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_chr.txt
sed -i 's/chr//g' z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_auto.txt
sed -i 's/chr//g' z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_sex.txt

#grab 1:1 ortholog transcripts from gff comparison file
grep 'transcript' gila2_galgal5_gff_comparison.txt > gila2_galgal5_gff_transcript.txt

#next command just excises two transcripts from gila2.refbased.stringtie_compiled_transcripts_sex.txt file, as they are in the PAR.
awk ' $5 !~ "304" || $6 >= 1700000 ' gila2_galgal5_gff_transcript.txt > gila2_galgal5_gff_transcripts.txt
sed -i 's/chr//g' gila2_galgal5_gff_transcripts.txt

#match Heloderma transcripts with 1:1 orthologs in Gallus and pull their associated stringtie info in Gallus
head -1 z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt > galgal5_gila2_stringtie_transcripts.txt
awk 'NR==FNR{c[$2$3]++;next};c[$2$3] > 0' gila2_galgal5_gff_transcripts.txt z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_chr.txt >> galgal5_gila2_stringtie_transcripts.txt

#extract transcripts from the representative autosomes and transcripts that aren't expressed in either sex
head -1 gila2.refbased.stringtie_compiled_per_transcript.txt > gila2.refbased.stringtie_compiled_transcripts_auto.txt
awk ' $2 == "0" || $2 == "1" || $2 == "2" || $2 == "3" ' gila2.refbased.stringtie_compiled_per_transcript.txt | grep -v "[[:space:]]0.0[[:space:]]" >> gila2.refbased.stringtie_compiled_transcripts_auto.txt

#extract transcripts from the putative sex-linked scaffolds and transcripts that aren't expressed in either sex
head -1 gila2.refbased.stringtie_compiled_per_transcript.txt > gila2.refbased.stringtie_compiled_transcripts_sex.txt
awk ' $2 == "157" || $2 == "218" || $2 == "304" || $2 == "398" ' gila2.refbased.stringtie_compiled_per_transcript.txt | grep -v "[[:space:]]0.0[[:space:]]" >> gila2.refbased.stringtie_compiled_transcripts_sex.txt

#match Heloderma transcripts with 1:1 orthologs in Gallus and pull their associated stringtie info in Heloderma
head -1 gila2.refbased.stringtie_compiled_per_transcript.txt > gila2_galgal5_stringtie_transcripts.txt
awk 'NR==FNR{c[$5$6]++;next};c[$2$3] > 0' gila2_galgal5_gff_transcripts.txt gila2.refbased.stringtie_compiled_transcripts_sex.txt >> gila2_galgal5_stringtie_transcripts.txt



#extract Heloderma genes from 1:1 list with Gallus, grouped by Gallus chromosome and match with gene expression
cat chr.txt | while read line 
do
head -1 gila_chicken_gff_comparison_FULL.txt | cut -f 1,5-6 > Heloderma_genes_Gallus_$line\.txt
awk -v awkvar="$line" ' $2 == 'awkvar' ' gila_chicken_gffcomparison_FULL.txt | cut -f 1,5-6 | uniq >> Heloderma_genes_Gallus_$line\.txt

head -1 gila2.refbased.stringtie_compiled_per_transcript.txt > gila2.refbased.stringtie_$line\.txt
awk 'NR==FNR{c[$2$3]++;next};c[$2$3] > 0' Heloderma_genes_Gallus_$line\.txt gila2.refbased.stringtie_compiled_per_transcript.txt | grep -v "[[:space:]]0.0[[:space:]]" >> gila2.refbased.stringtie_$line\.txt
done

head -1 z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt > z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_all.txt
grep "chr" z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt | grep -v "_" | grep -v "[[:space:]]0.0[[:space:]]" | grep -v "M" | grep -v "W" | grep -v "LGE64" | grep -v "chr16">> z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_all.txt

