import os

configfile: "gilas_config.json"

fastq_directory = "/mnt/storage/SAYRES/GilaMonster/all_fastqs/"
temp_directory = "temp/"

bbduksh_path = "bbduk.sh"
bwa_path = "bwa"
fastqc_path = "fastqc"
multiqc_path = "multiqc"
samblaster_path = "samblaster"
samtools_path = "samtools"

samples = [
	"G_10_dna", "G_10_rna", "G_16_dna", "G_16_rna", "G_30_dna", "G_30_rna",
	"G_35_dna", "G_35_rna", "G_KI01_dna", "G_KI01_rna", "G_L_dna", "G_L_rna"]

dna = ["G_10_dna", "G_16_dna", "G_30_dna", "G_35_dna", "G_KI01_dna", "G_L_dna"]

rna = ["G_10_rna", "G_16_rna", "G_30_rna", "G_35_rna", "G_KI01_rna", "G_L_rna"]

fastq_prefixes = [config[x]["fq1"][:-9] for x in samples] + [config[x]["fq2"][:-9] for x in samples]

rule all:
	input:
		expand(
			"new_reference/{assembly}.fasta.fai",
			assembly=["gila1"]),
		"multiqc/multiqc_report.html",
		"multiqc_trimmed_dna/multiqc_report.html",
		expand(
			"processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai",
			sample=dna, genome=["gila1"]),
		expand(
			"stats/{sample}.{genome}.dna.mkdup.sorted.bam.stats",
			sample=dna, genome=["gila1"]),

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
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o multiqc fastqc"

rule trim_adapters_paired_bbduk_dna:
	input:
		fq1 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq1"]),
		fq2 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq2"])
	output:
		out_fq1 = "trimmed_dna_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_dna_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=75 maq=20"

rule fastqc_analysis_trimmed_dna:
	input:
		"trimmed_dna_fastqs/{sample}_trimmed_{read}.fastq.gz"
	output:
		"fastqc_trimmed_dna/{sample}_trimmed_{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed_dna {input}"

rule multiqc_analysis_trimmed_dna:
	input:
		expand(
			"fastqc_trimmed_dna/{sample}_trimmed_{read}_fastqc.html",
			sample=dna, read=["read1", "read2"])
	output:
		"multiqc_trimmed_dna/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o multiqc_trimmed_dna fastqc_trimmed_dna"

rule map_and_process_trimmed_reads:
	input:
		fq1 = "trimmed_dna_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_dna_fastqs/{sample}_trimmed_read2.fastq.gz",
		amb = "new_reference/{genome}.fasta.amb",
		ref = "new_reference/{genome}.fasta"
	output:
		"processed_bams/{sample}.{genome}.mkdup.sorted.bam"
	params:
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"],
		bwa = bwa_path,
		samblaster = samblaster_path,
		samtools = samtools_path,
		threads = 4
	threads: 4
	shell:
		"{params.bwa} mem -t {params.threads} -R "
		"'@RG\\tID:{params.id}\\tSM:{params.sm}\\tLB:{params.lb}\\tPU:{params.pu}\\tPL:{params.pl}' "
		"{input.ref} {input.fq1} {input.fq2}"
		"| {params.samblaster} | {params.samtools} fixmate -O bam - - "
		"| {params.samtools} sort -O bam -o {output}"

rule index_bam:
	input:
		"processed_bams/{sample}.{genome}.mkdup.sorted.bam"
	output:
		"processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} index {input}"

rule bam_stats_dna:
	input:
		bam = "processed_bams/{sample}.{genome}.mkdup.sorted.bam",
		bai = "processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai"
	output:
		"stats/{sample}.{genome}.dna.mkdup.sorted.bam.stats"
	params:
		samtools = samtools_path
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"
