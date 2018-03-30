import os

configfile: "gilas_config.json"

fastq_directory = "/mnt/storage/SAYRES/GilaMonster/all_fastqs/"
temp_directory = "temp/"

bbduksh_path = "bbduk.sh"
bwa_path = "bwa"
fastqc_path = "fastqc"
gatk_path = "gatk-launch"
multiqc_path = "multiqc"
samblaster_path = "samblaster"
samtools_path = "samtools"
xyalign_path = "xyalign"
hisat2_build_path = "hisat2-build"
xyalign_anaconda_env = "xyalign_env"

samples = [
	"G_10_dna", "G_10_rna", "G_16_dna", "G_16_rna", "G_30_dna", "G_30_rna",
	"G_35_dna", "G_35_rna", "G_KI01_dna", "G_KI01_rna", "G_L_dna", "G_L_rna"]

dna = ["G_10_dna", "G_16_dna", "G_30_dna", "G_35_dna", "G_KI01_dna", "G_L_dna"]

rna = ["G_10_rna", "G_16_rna", "G_30_rna", "G_35_rna", "G_KI01_rna", "G_L_rna"]

fastq_prefixes = [
	config[x]["fq1"][:-9] for x in samples] + [
		config[x]["fq2"][:-9] for x in samples]

# Parameters for splitting reference
num_chunks = 25
chunk_range = [x for x in range(1, num_chunks + 1)]

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
		expand("hisat2_index/{assembly}.8.ht2", assembly=["gila1"]),
		"multiqc_trimmed_rna/multiqc_report.html",
		expand(
			"xyalign_analyses/{genome}/results/{genome}_chrom_stats_count.txt",
			genome=["gila1"]),
		# expand(
		# 	"vcf/{sample}.{genome}.{chunk}.g.vcf.gz",
		# 	sample=dna, genome=["gila1"], chunk=chunk_range),
		expand(
			"vcf/{genome}.{chunk}.gatk.raw.vcf.gz",
			genome=["gila1"], chunk=chunk_range)

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

rule chunk_reference:
	input:
		fai = "new_reference/{assembly}.fasta.fai"
	output:
		expand("new_reference/{{assembly}}_split_chunk{num}.bed", num=chunk_range)
	params:
		chunks = num_chunks
	shell:
		"python scripts/Chunk_fai.py --fai {input.fai} "
		"--out_prefix new_reference/split --chunks {params.chunks}"

rule prepare_hisat_index:
	input:
		"new_reference/{assembly}.fasta"
	output:
		"hisat2_index/{assembly}.8.ht2"
	params:
		hisat2_build = hisat2_build_path,
		prefix = "hisat2_index/{assembly}"
	shell:
		"{params.hisat2_build} {input} {params.prefix}"

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

rule trim_adapters_paired_bbduk_rna:
	input:
		fq1 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq1"]),
		fq2 = lambda wildcards: os.path.join(
			fastq_directory, config[wildcards.sample]["fq2"])
	output:
		out_fq1 = "trimmed_rna_fastqs/{sample}_trimmed_read1.fastq.gz",
		out_fq2 = "trimmed_rna_fastqs/{sample}_trimmed_read2.fastq.gz"
	params:
		bbduksh = bbduksh_path
	shell:
		"{params.bbduksh} -Xmx3g in1={input.fq1} in2={input.fq2} "
		"out1={output.out_fq1} out2={output.out_fq2} "
		"ref=misc/adapter_sequence.fa ktrim=r k=21 mink=11 hdist=2 tbo tpe "
		"qtrim=rl trimq=15 minlen=60 maq=20"

rule fastqc_analysis_trimmed_rna:
	input:
		"trimmed_rna_fastqs/{sample}_trimmed_{read}.fastq.gz"
	output:
		"fastqc_trimmed_rna/{sample}_trimmed_{read}_fastqc.html"
	params:
		fastqc = fastqc_path
	shell:
		"{params.fastqc} -o fastqc_trimmed_rna {input}"

rule multiqc_analysis_trimmed_rna:
	input:
		expand(
			"fastqc_trimmed_rna/{sample}_trimmed_{read}_fastqc.html",
			sample=rna, read=["read1", "read2"])
	output:
		"multiqc_trimmed_rna/multiqc_report.html"
	params:
		multiqc = multiqc_path
	shell:
		"export LC_ALL=en_US.UTF-8 && export LANG=en_US.UTF-8 && "
		"{params.multiqc} -o multiqc_trimmed_rna fastqc_trimmed_rna"

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

rule chrom_stats_dna:
	input:
		bams = lambda wildcards: expand(
			"processed_bams/{sample}.{genome}.mkdup.sorted.bam",
			sample=dna, genome=[wildcards.genome]),
		bais = lambda wildcards: expand(
			"processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai",
			sample=dna, genome=[wildcards.genome])
	output:
		counts = "xyalign_analyses/{genome}/results/{genome}_chrom_stats_count.txt"
	params:
		xyalign = xyalign_path,
		sample_id = "{genome}",
		xyalign_env = xyalign_anaconda_env
	shell:
		"source activate {params.xyalign_env} && "
		"{params.xyalign} --CHROM_STATS --use_counts "
		"--chromosomes ALL --bam {input.bams} --ref null "
		"--sample_id {params.sample_id} "
		"--output_dir xyalign_analyses/{params.sample_id}"

rule gatk_gvcf_per_chunk:
	input:
		ref = "new_reference/{genome}.fasta",
		bam = "processed_bams/{sample}.{genome}.mkdup.sorted.bam",
		bai = "processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai",
		chunkfile = "new_reference/{genome}_split_chunk{chunk}.bed"
	output:
		"vcf/{sample}.{genome}.{chunk}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path
	threads:
		4
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.chunkfile} """
		"""-ERC GVCF -O {output}"""

rule gatk_combinegvcfs_per_chunk:
	input:
		ref = "new_reference/{genome}.fasta",
		gvcfs = expand(
			"vcf/{sample}.{{genome}}.{{chunk}}.g.vcf.gz", sample=dna)
	output:
		"vcf/{genome}.{chunk}.gatk.raw.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path
	threads:
		4
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}""")
