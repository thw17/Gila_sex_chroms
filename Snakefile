import os

configfile: "gilas_config.json"

fastq_directory = "/mnt/storage/SAYRES/GilaMonster/all_fastqs/"
temp_directory = "temp/"

bwa_path = "bwa"
fastqc_path = "fastqc"
multiqc_path = "multiqc"
samtools_path = "samtools"

samples = [
	"G_10_dna", "G_10_rna", "G_16_dna", "G_16_rna", "G_30_dna", "G_30_rna",
	"G_35_dna", "G_35_rna", "G_KI01_dna", "G_KI01_rna", "G_L_dna", "G_L_rna"]

fastq_prefixes = [config[x]["fq1"][:-9] for x in samples] + [config[x]["fq2"][:-9] for x in samples]

rule all:
	input:
		expand("new_reference/{assembly}.fasta.fai", assembly=["gila1"]),
		"multiqc/multiqc_report.html"


rule prepare_reference:
	input:
		ref = lambda wildcards: config["genome_paths"][wildcards.assembly]
	output:
		new = "new_reference/{assembly}.fasta",
		fai = "new_reference/{assembly}.fasta.fai",
		amb = "new_reference/{assembly}.fasta.amb",
		dict = "new_reference/{assembly}.dict"
	params:
		samtools = samtools_path,
		bwa = bwa_path
	run:
		shell(
			"ln -s ../{} {{output.new}} && touch -h {{output.new}}".format(input.ref))
		# faidx
		shell(
			"{params.samtools} faidx {output.new}")
		# .dict
		shell(
			"{params.samtools} dict -o {output.dict} {output.new}")
		# bwa
		shell(
			"{params.bwa} index {output.new}")

rule fastqc_analysis:
	input:
		os.path.join(fastq_directory, "{fq_prefix}.fastq.gz")
	output:
		"fastqc/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc {input}"

rule multiqc_analysis:
	input:
		expand("fastqc/{fq_prefix}_fastqc.html", fq_prefix=fastq_prefixes)
	output:
		"multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"{params.multiqc} -o multiqc fastqc"
