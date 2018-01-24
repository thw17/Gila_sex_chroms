import os

configfile: "gilas_config.json"

fastq_directory = "/mnt/storage/SAYRES/GilaMonster/all_fastqs/"
temp_directory = "temp/"

bwa_path = "bwa"
fastqc_path = "fastqc"
samtools_path = "samtools"

rule all:
	input:
		expand("reference/{assembly}.fasta.fai", assembly=["gila1"])


rule prepare_reference:
	input:
		ref = lambda wildcards: config["genome_paths"][wildcards.assembly]
	output:
		new = "reference/{assembly}.fasta",
		fai = "reference/{assembly}.fasta.fai",
		amb = "reference/{assembly}.fasta.amb",
		dict = "reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		shell("ln -s ../{} {{output.new}} && touch -h {{output.new}}".format(input.ref))
		# faidx
		shell("{params.samtools} faidx {input}")
		# .dict
		shell("{params.samtools} dict -o {output.dict} {input}")
		# bwa
		shell("{params.bwa} index {input}")

# rule fastqc_analysis:
# 	input:
#
# 	output:
# 		"fastqc/{fq_prefix}_fastqc.html"
# 	params:
# 		fastqc = fastqc_path
# 	shell:
# 		"{params.fastqc} -o fastqc {input}"
