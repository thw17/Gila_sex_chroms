#!bin/bash
#3 files needed in the current directory for this code to run:
#"z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt.gz"
#"gila2_galgal5_gff_comparison.txt.gz"
#"gila2.refbased.stringtie_compiled_per_transcript.txt.gz"

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

#next command just excises two transcripts from gila2.refbased.stringtie_compiled_transcripts_sex.txt file as they are in the PAR.
#10688	304	366570	4.982189	5.659705666666667	1.1359877488924381
#10689	304	752916	124.99868299999999	153.97120166666664	1.2317825913947162

awk ' $5 !~ "304" || $6 >= 1700000 ' gila2_galgal5_gff_transcript.txt > gila2_galgal5_gff_transcripts.txt
sed -i 's/chr//g' gila2_galgal5_gff_transcripts.txt

#match Heloderma transcripts with 1:1 orthologs in Gallus and pull their associated stringtie info in Gallus
head -1 z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript.txt > galgal5_gila2_stringtie_transcripts.txt
awk 'NR==FNR{c[$2$3]++;next};c[$2$3] > 0' gila2_galgal5_gff_transcripts.txt z_ortho_filtered-gila2_galgal5-galgal5.refbased.stringtie_compiled_per_transcript_chr.txt >> galgal5_gila2_stringtie_transcripts.txt

#extract transcripts from the representative autosomes and transcripts that aren't expressed in either sex
head -1 gila2.refbased.stringtie_compiled_per_transcript.txt > gila2.refbased.stringtie_compiled_transcripts_auto.txt
awk ' $2 == "0" ' gila2.refbased.stringtie_compiled_per_transcript.txt | grep -v "[[:space:]]0.0[[:space:]]" >> gila2.refbased.stringtie_compiled_transcripts_auto.txt
awk ' $2 == "1" ' gila2.refbased.stringtie_compiled_per_transcript.txt | grep -v "[[:space:]]0.0[[:space:]]" >> gila2.refbased.stringtie_compiled_transcripts_auto.txt
awk ' $2 == "2" ' gila2.refbased.stringtie_compiled_per_transcript.txt | grep -v "[[:space:]]0.0[[:space:]]" >> gila2.refbased.stringtie_compiled_transcripts_auto.txt
awk ' $2 == "3" ' gila2.refbased.stringtie_compiled_per_transcript.txt | grep -v "[[:space:]]0.0[[:space:]]" >> gila2.refbased.stringtie_compiled_transcripts_auto.txt

#extract transcripts from the putative sex-linked scaffolds and transcripts that aren't expressed in either sex
head -1 gila2.refbased.stringtie_compiled_per_transcript.txt > gila2.refbased.stringtie_compiled_transcripts_sex.txt
awk ' $2 == "157" ' gila2.refbased.stringtie_compiled_per_transcript.txt | grep -v "[[:space:]]0.0[[:space:]]" >> gila2.refbased.stringtie_compiled_transcripts_sex.txt
awk ' $2 == "218" ' gila2.refbased.stringtie_compiled_per_transcript.txt | grep -v "[[:space:]]0.0[[:space:]]" >> gila2.refbased.stringtie_compiled_transcripts_sex.txt
awk ' $2 == "304" ' gila2.refbased.stringtie_compiled_per_transcript.txt | grep -v "[[:space:]]0.0[[:space:]]" >> gila2.refbased.stringtie_compiled_transcripts_sex.txt
awk ' $2 == "398" ' gila2.refbased.stringtie_compiled_per_transcript.txt | grep -v "[[:space:]]0.0[[:space:]]" >> gila2.refbased.stringtie_compiled_transcripts_sex.txt

#match Heloderma transcripts with 1:1 orthologs in Gallus and pull their associated stringtie info in Heloderma
head -1 gila2.refbased.stringtie_compiled_per_transcript.txt > gila2_galgal5_stringtie_transcripts.txt
awk 'NR==FNR{c[$5$6]++;next};c[$2$3] > 0' gila2_galgal5_gff_transcripts.txt gila2.refbased.stringtie_compiled_transcripts_sex.txt >> gila2_galgal5_stringtie_transcripts.txt









#additional files needed for 'corrected' re-analysis
#"z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased.txt"
#"corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals.txt"


#replicate for "corrected" data
#cut -f 1-3,13-15 z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased.txt > z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset.txt
#head -1 z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset.txt > z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_chr.txt
#head -1 z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset.txt > z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_auto.txt
#head -1 z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset.txt > z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_sex.txt
#grep "chr"  z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset.txt | grep -v "_" | grep -v "[[:space:]]0.0[[:space:]]" | grep -v "M" | grep -v "W" | grep -v "LGE64" >> z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_chr.txt
#grep "chr"  z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset.txt | grep -v "_" | grep -v "[[:space:]]0.0[[:space:]]" | grep -v "28" | grep -v "M" | grep -v "W" | grep -v "Z" | grep -v "LGE64" >> z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_auto.txt
#grep "chrZ" z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset.txt | grep -v "_" | grep -v "[[:space:]]0.0[[:space:]]" >> z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_sex.txt
#sed -i 's/chr//g' z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_chr.txt
#sed -i 's/chr//g' z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_auto.txt
#sed -i 's/chr//g' z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_sex.txt
#
##match Heloderma transcripts with 1:1 orthologs in Gallus and pull their associated stringtie info in Gallus
#head -1 z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_chr.txt > galgal5_gila2_stringtie_transcripts_corrected.txt
#awk 'NR==FNR{c[$2$3]++;next};c[$2$3] > 0' gila2_galgal5_gff_transcripts.txt z_ortho_filtered.corrected.stringtie_compiled_per_transcript_separate_individuals.gila2_galgal5.galgal5_refbased_subset_chr.txt >> galgal5_gila2_stringtie_transcripts_corrected.txt
#
#cut -f 1-3,13-15 corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals.txt > corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset.txt
#head -1 corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset.txt > corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_auto.txt
#awk ' $2 == "0" ' corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset.txt | grep -v "[[:space:]]0.0[[:space:]]" >> corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_auto.txt
#awk ' $2 == "1" ' corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset.txt | grep -v "[[:space:]]0.0[[:space:]]" >> corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_auto.txt
#awk ' $2 == "2" ' corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset.txt | grep -v "[[:space:]]0.0[[:space:]]" >> corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_auto.txt
#awk ' $2 == "3" ' corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset.txt | grep -v "[[:space:]]0.0[[:space:]]" >> corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_auto.txt
#
#head -1 corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset.txt > corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_sex.txt
#awk ' $2 == "157" ' corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset.txt | grep -v "[[:space:]]0.0[[:space:]]" >> corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_sex.txt
#awk ' $2 == "218" ' corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset.txt | grep -v "[[:space:]]0.0[[:space:]]" >> corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_sex.txt
#awk ' $2 == "304" ' corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset.txt | grep -v "[[:space:]]0.0[[:space:]]" >> corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_sex.txt
#awk ' $2 == "398" ' corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset.txt | grep -v "[[:space:]]0.0[[:space:]]" >> corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_sex.txt
#
##Excise these two rows from gila2.refbased.stringtie_compiled_transcripts_sex.txt file as they are in the PAR.
##10688	304	366570	4.982189	5.659705666666667	1.1359877488924381
##10689	304	752916	124.99868299999999	153.97120166666664	1.2317825913947162
#
##match Heloderma transcripts with 1:1 orthologs in Gallus and pull their associated stringtie info in Heloderma
#head -1 corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_sex.txt > gila2_galgal5_stringtie_transcripts_corrected.txt
#awk 'NR==FNR{c[$5$6]++;next};c[$2$3] > 0' gila2_galgal5_gff_transcripts.txt corrected.gila2.refbased.stringtie_compiled_per_transcript_separate_individuals_subset_sex.txt >> gila2_galgal5_stringtie_transcripts_corrected.txt
#








