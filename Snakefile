import os

configfile: "gilas_config.json"

# fastq_directory = "/mnt/storage/SAYRES/GilaMonster/all_fastqs/"
fastq_directory = "fastqs"
temp_directory = "temp/"

# Runtime values
very_short = "6:00:00"
medium = "12:00:00"
day = "24:00:00"
long = "48:00:00"
very_long = "72:00:00"

bbduksh_path = "bbduk.sh"
bcftools_path = "bcftools"
bwa_path = "bwa"
fastqc_path = "fastqc"
gatk_path = "gatk"
bgzip_path = "bgzip"
multiqc_path = "multiqc"
samblaster_path = "samblaster"
samtools_path = "samtools"
tabix_path = "tabix"
xyalign_path = "xyalign"
hisat2_path = "hisat2"
hisat2_build_path = "hisat2-build"
stringtie_path = "stringtie"
xyalign_anaconda_env = "gila_xyalign_env"

lastal_path = "lastal"
lastdb_path = "lastdb"

# assembly_list = ["gila1", "gila2"]
assembly_list = ["gila2"]
scaffolds_to_analyze = ["157", "218", "304", "398", "674", "0", "1", "2", "3"]

samples = [
	"G_10_dna", "G_10_rna", "G_16_dna", "G_16_rna", "G_30_dna", "G_30_rna",
	"G_35_dna", "G_35_rna", "G_KI01_dna", "G_KI01_rna", "G_L_dna", "G_L_rna"]

dna = ["G_10_dna", "G_16_dna", "G_30_dna", "G_35_dna", "G_KI01_dna", "G_L_dna"]

rna = ["G_10_rna", "G_16_rna", "G_30_rna", "G_35_rna", "G_KI01_rna", "G_L_rna"]

rna_to_dna = {
	"G_10_rna": "G_10_dna",
	"G_16_rna": "G_16_dna",
	"G_30_rna": "G_30_dna",
	"G_35_rna": "G_35_dna",
	"G_KI01_rna": "G_KI01_dna",
	"G_L_rna": "G_L_dna"
}

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
			assembly=assembly_list),
		"multiqc/multiqc_report.html",
		"multiqc_trimmed_dna/multiqc_report.html",
		expand(
			"processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai",
			sample=dna, genome=assembly_list),
		expand(
			"stats/{sample}.{genome}.dna.mkdup.sorted.bam.stats",
			sample=dna, genome=assembly_list),
		expand(
			"xyalign_analyses/{genome}/results/{genome}_chrom_stats_count.txt",
			genome=assembly_list),
		# expand(
		# 	"vcf/{sample}.{genome}.{chunk}.g.vcf.gz",
		# 	sample=dna, genome=["gila1"], chunk=chunk_range),
		# expand(
		# 	"combined_gvcfs/{genome}.{chunk}.gatk.combinegvcf.g.vcf.gz",
		# 	genome=["gila1"], chunk=chunk_range),
		expand(
		 	"xyalign_analyses/{sample}.{genome}/logfiles/{sample}.{genome}_xyalign.log",
		 	sample=dna, genome=assembly_list),
		expand(
			"genotyped_vcfs/{genome}.{chunk}.gatk.called.raw.vcf.gz",
			genome=assembly_list, chunk=chunk_range),
		expand(
			"stringtie_gtfs_{strategy}/{sample}_{genome}/{sample}.{genome}.secondpass.gtf",
			strategy=["mixed", "denovo", "refbased"], genome=assembly_list, sample=rna),
		expand(
			"stats/{sample}.{genome}.rna.sorted.bam.stats",
			genome=assembly_list, sample=rna),
		expand(
			"results/{genome}.{strategy}.stringtie_compiled_per_transcript.txt",
		 	strategy=["mixed", "denovo", "refbased"],
		 	genome=assembly_list),
		# expand(
		# 	"results/{genome}.chromstats_compiled.txt",
		# 	genome=assembly_list),
		# expand(
		# 	"combined_vcfs/combined.{genome}.filtered.vcf.gz.tbi",
		# 	genome=assembly_list),
		# expand(
		# 	"results/{genome}.het_rate.txt",
		# 	genome=assembly_list),
		expand(
			"results/all_compiled.{genome}.{strategy}.txt",
			strategy=["mixed", "denovo", "refbased"],
			genome=assembly_list),
		expand(
			"komodo/komodo_scaff218_{assembly}_align.maf",
			assembly=assembly_list),
		expand(
			"par_results/scaffold{scaff}_{genome}.txt",
			genome=assembly_list,
			scaff=scaffolds_to_analyze)

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
		bwa = bwa_path,
		threads = 4,
		mem = 16,
		t = medium
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
		chunks = num_chunks,
		out_prefix = "new_reference/{assembly}_split",
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"python scripts/Chunk_fai.py --fai {input.fai} "
		"--out_prefix {params.out_prefix} --chunks {params.chunks}"

rule hisat2_reference_index:
	input:
		"new_reference/{genome}.fasta"
	output:
		expand(
			"new_reference/hisat2/{{genome}}.{suffix}.ht2",
			suffix=[
				"1", "2", "3", "4", "5", "6", "7", "8"])
	params:
		hisat2_build = hisat2_build_path,
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"{params.hisat2_build} {input} new_reference/hisat2/{wildcards.genome}"

rule hisat2_map_reads:
	input:
		idx = expand(
			"new_reference/hisat2/{{genome}}.{suffix}.ht2",
			suffix=["1", "2", "3", "4", "5", "6", "7", "8"]),
		fq1 = "trimmed_rna_fastqs/{sample}_trimmed_read1.fastq.gz",
		fq2 = "trimmed_rna_fastqs/{sample}_trimmed_read2.fastq.gz"
	output:
		"processed_rna_bams/{sample}.{genome}.sorted.bam"
	params:
		threads = 4,
		mem = 16,
		t = long,
		hisat2 = hisat2_path,
		samtools = samtools_path,
		id = lambda wildcards: config[wildcards.sample]["ID"],
		sm = lambda wildcards: config[wildcards.sample]["SM"],
		lb = lambda wildcards: config[wildcards.sample]["LB"],
		pu = lambda wildcards: config[wildcards.sample]["PU"],
		pl = lambda wildcards: config[wildcards.sample]["PL"]
	shell:
		"{params.hisat2} -p {params.threads} --dta "
		"--rg-id {params.id} --rg SM:{params.sm} --rg LB:{params.lb} "
		"--rg PU:{params.pu} --rg PL:{params.pl} "
		"-x new_reference/hisat2/{wildcards.genome} "
		"-1 {input.fq1} -2 {input.fq2} | "
		"{params.samtools} sort -O bam -o {output}"

rule index_bam_rna:
	input:
		"processed_rna_bams/{sample}.{genome}.sorted.bam"
	output:
		"processed_rna_bams/{sample}.{genome}.sorted.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule bam_stats_rna:
	input:
		bam = "processed_rna_bams/{sample}.{genome}.sorted.bam",
		bai = "processed_rna_bams/{sample}.{genome}.sorted.bam.bai"
	output:
		"stats/{sample}.{genome}.rna.sorted.bam.stats"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} stats {input.bam} | grep ^SN | cut -f 2- > {output}"

rule stringtie_first_pass_denovo:
	input:
		bam = "processed_rna_bams/{sample}.{genome}.sorted.bam"
	output:
		"stringtie_gtfs_denovo/{sample}_{genome}/{sample}.{genome}.firstpass.gtf"
	threads:
		4
	params:
		stringtie = stringtie_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.stringtie} {input.bam} -o {output} -p {params.threads}"

rule create_stringtie_merged_list_denovo:
	input:
		lambda wildcards: expand(
			"stringtie_gtfs_denovo/{sample}_{genome}/{sample}.{genome}.firstpass.gtf",
			genome=wildcards.genome,
			sample=rna)
	output:
		"stringtie_gtfs_denovo/{genome}_gtflist.txt"
	params:
		threads = 1,
		mem = 4,
		t = very_short
	run:
		shell("echo -n > {output}")
		for i in input:
			shell("echo {} >> {{output}}".format(i))

rule stringtie_merge_denovo:
	input:
		bam_list = "stringtie_gtfs_denovo/{genome}_gtflist.txt"
	output:
		"stringtie_gtfs_denovo/{genome}.merged.gtf"
	threads:
		4
	params:
		stringtie = stringtie_path,
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"{params.stringtie} --merge {input.bam_list} -o {output} -p {params.threads}"

rule stringtie_second_pass_denovo:
	input:
		bam = "processed_rna_bams/{sample}.{genome}.sorted.bam",
		gtf = "stringtie_gtfs_denovo/{genome}.merged.gtf"
	output:
		gtf = "stringtie_gtfs_denovo/{sample}_{genome}/{sample}.{genome}.secondpass.gtf",
		ctab = "stringtie_gtfs_denovo/{sample}_{genome}/t_data.ctab"
	threads:
		4
	params:
		stringtie = stringtie_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.stringtie} {input.bam} -o {output.gtf} -p {params.threads} "
		"-G {input.gtf} -B -e"

rule stringtie_refbased:
	input:
		bam = "processed_rna_bams/{sample}.{genome}.sorted.bam",
		gff = lambda wildcards: config["annotation"][wildcards.genome]
	output:
		gtf = "stringtie_gtfs_refbased/{sample}_{genome}/{sample}.{genome}.secondpass.gtf",
		ctab = "stringtie_gtfs_refbased/{sample}_{genome}/t_data.ctab"
	threads:
		4
	params:
		stringtie = stringtie_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.stringtie} {input.bam} -o {output.gtf} -p {params.threads} "
		"-G {input.gff} -B -e"

rule stringtie_first_pass_mixed:
	input:
		bam = "processed_rna_bams/{sample}.{genome}.sorted.bam",
		gff = lambda wildcards: config["annotation"][wildcards.genome]
	output:
		"stringtie_gtfs_mixed/{sample}_{genome}/{sample}.{genome}.firstpass.gtf"
	params:
		stringtie = stringtie_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.stringtie} {input.bam} -o {output} -p {params.threads} "
		"-G {input.gff}"

rule create_stringtie_merged_list_mixed:
	input:
		lambda wildcards: expand(
			"stringtie_gtfs_mixed/{sample}_{genome}/{sample}.{genome}.firstpass.gtf",
			genome=wildcards.genome,
			sample=rna)
	output:
		"stringtie_gtfs_mixed/{genome}_gtflist.txt"
	params:
		threads = 1,
		mem = 4,
		t = very_short
	run:
		shell("echo -n > {output}")
		for i in input:
			shell("echo {} >> {{output}}".format(i))

rule stringtie_merge_mixed:
	input:
		bam_list = "stringtie_gtfs_mixed/{genome}_gtflist.txt"
	output:
		"stringtie_gtfs_mixed/{genome}.merged.gtf"
	params:
		stringtie = stringtie_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.stringtie} --merge {input.bam_list} -o {output} -p {params.threads}"

rule stringtie_second_pass_mixed:
	input:
		bam = "processed_rna_bams/{sample}.{genome}.sorted.bam",
		gtf = "stringtie_gtfs_mixed/{genome}.merged.gtf"
	output:
		gtf = "stringtie_gtfs_mixed/{sample}_{genome}/{sample}.{genome}.secondpass.gtf",
		ctab = "stringtie_gtfs_mixed/{sample}_{genome}/t_data.ctab"
	params:
		stringtie = stringtie_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.stringtie} {input.bam} -o {output.gtf} -p {params.threads} "
		"-G {input.gtf} -B -e"

rule fastqc_analysis:
	input:
		os.path.join(fastq_directory, "{fq_prefix}.fastq.gz")
	output:
		"fastqc/{fq_prefix}_fastqc.html"
	params:
		fastqc = fastqc_path,
		threads = 1,
		mem = 4,
		t = very_short
	shell:
		"{params.fastqc} -o fastqc {input}"

rule multiqc_analysis:
	input:
		expand("fastqc/{fq_prefix}_fastqc.html", fq_prefix=fastq_prefixes)
	output:
		"multiqc/multiqc_report.html"
	params:
		multiqc = multiqc_path,
		threads = 1,
		mem = 4,
		t = very_short
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
		bbduksh = bbduksh_path,
		threads = 2,
		mem = 8,
		t = very_short
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
		fastqc = fastqc_path,
		threads = 1,
		mem = 4,
		t = very_short
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
		multiqc = multiqc_path,
		threads = 1,
		mem = 4,
		t = very_short
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
		bbduksh = bbduksh_path,
		threads = 2,
		mem = 8,
		t = very_short
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
		fastqc = fastqc_path,
		threads = 1,
		mem = 4,
		t = very_short
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
		multiqc = multiqc_path,
		threads = 1,
		mem = 4,
		t = very_short
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
		threads = 4,
		mem = 16,
		t = very_long
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
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule bam_stats_dna:
	input:
		bam = "processed_bams/{sample}.{genome}.mkdup.sorted.bam",
		bai = "processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai"
	output:
		"stats/{sample}.{genome}.dna.mkdup.sorted.bam.stats"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
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
		xyalign_env = xyalign_anaconda_env,
		threads = 4,
		mem = 16,
		t = long
	conda:
		"envs/gila_xyalign_environment.yaml"
	shell:
		"{params.xyalign} --CHROM_STATS --use_counts "
		"--chromosomes ALL --bam {input.bams} --ref null "
		"--sample_id {params.sample_id} "
		"--output_dir xyalign_analyses/{params.sample_id}"

rule bam_analysis_dna:
	input:
		bam = "processed_bams/{sample}.{genome}.mkdup.sorted.bam",
		bai = "processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai",
		ref = "new_reference/{genome}.fasta",
		fai = "new_reference/{genome}.fasta.fai",
	output:
		bed = "xyalign_analyses/{sample}.{genome}/bed/{sample}.{genome}_full_dataframe_depth_mapq_preprocessing.csv",
		log = "xyalign_analyses/{sample}.{genome}/logfiles/{sample}.{genome}_xyalign.log"
	params:
		xyalign = xyalign_path,
		sample_id = "{sample}.{genome}",
		xyalign_env = xyalign_anaconda_env,
		threads = 4,
		mem = 16,
		t = long
	conda:
		"envs/gila_xyalign_environment.yaml"
	shell:
		"{params.xyalign} --ANALYZE_BAM "
		"--chromosomes 0 1 2 3 157 218 304 398 "
		"--bam {input.bam} --ref {input.ref} "
		"--sample_id {params.sample_id} "
		"--output_dir xyalign_analyses/{params.sample_id} "
		"--cpus 4 --window_size 5000"

rule gatk_gvcf_per_chunk:
	input:
		ref = "new_reference/{genome}.fasta",
		bam = "processed_bams/{sample}.{genome}.mkdup.sorted.bam",
		bai = "processed_bams/{sample}.{genome}.mkdup.sorted.bam.bai",
		chunkfile = "new_reference/{genome}_split_chunk{chunk}.bed"
	output:
		"gvcfs/{sample}.{genome}.{chunk}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 4,
		mem = 16,
		t = very_long
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} -L {input.chunkfile} """
		"""-ERC GVCF --do-not-run-physical-phasing -O {output}"""

rule gatk_combinegvcfs_per_chunk:
	input:
		ref = "new_reference/{genome}.fasta",
		gvcfs = expand(
			"gvcfs/{sample}.{{genome}}.{{chunk}}.g.vcf.gz", sample=dna)
	output:
		"combined_gvcfs/{genome}.{chunk}.gatk.combinegvcf.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 4,
		mem = 16,
		t = very_long
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}""")

rule gatk_genotypegvcf_per_chunk:
	input:
		ref = "new_reference/{genome}.fasta",
		gvcf = "combined_gvcfs/{genome}.{chunk}.gatk.combinegvcf.g.vcf.gz"
	output:
		"genotyped_vcfs/{genome}.{chunk}.gatk.called.raw.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 4,
		mem = 16,
		t = very_long
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output}"""

rule concatenate_split_vcfs:
	input:
		vcf = lambda wildcards: expand(
			"genotyped_vcfs/{gen}.{chunk}.gatk.called.raw.vcf.gz",
			gen=wildcards.genome,
			chunk=chunk_range)
	output:
		"combined_vcfs/combined.{genome}.raw.vcf.gz"
	params:
		bcftools = bcftools_path,
		threads = 2,
		mem = 8,
		t = medium
	shell:
		"{params.bcftools} concat -O z -o {output} {input.vcf}"

rule index_concatenated_vcf:
	input:
		"combined_vcfs/combined.{genome}.raw.vcf.gz"
	output:
		"combined_vcfs/combined.{genome}.raw.vcf.gz.tbi"
	params:
		tabix = tabix_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.tabix} -p vcf {input}"

rule filter_concatenated_vcf:
	input:
		vcf = "combined_vcfs/combined.{genome}.raw.vcf.gz",
		idx = "combined_vcfs/combined.{genome}.raw.vcf.gz.tbi"
	output:
		"combined_vcfs/combined.{genome}.filtered.vcf.gz"
	params:
		bgzip = bgzip_path,
		bcftools = bcftools_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.bcftools} filter -i "
		"'QUAL >= 30 && MQ >= 30 && QD > 2' {input.vcf} | "
		"{params.bcftools} filter -i 'FMT/DP >= 10 & FMT/GQ >= 30' -S . - | "
		"{params.bgzip} > {output}"

rule index_filtered_vcf:
	input:
		"combined_vcfs/combined.{genome}.filtered.vcf.gz"
	output:
		"combined_vcfs/combined.{genome}.filtered.vcf.gz.tbi"
	params:
		tabix = tabix_path,
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"{params.tabix} -p vcf {input}"

rule create_rna_bam_header:
	input:
		rbam = "processed_rna_bams/{sample}.{genome}.sorted.bam",
		rbai = "processed_rna_bams/{sample}.{genome}.sorted.bam.bai",
		gbam = lambda wildcards: expand(
			"processed_bams/{gsample}.{assembly}.mkdup.sorted.bam",
			gsample=[rna_to_dna[wildcards.sample]],
			assembly=wildcards.genome),
		gbai = lambda wildcards: expand(
			"processed_bams/{gsample}.{assembly}.mkdup.sorted.bam.bai",
			gsample=[rna_to_dna[wildcards.sample]],
			assembly=wildcards.genome)
	output:
		"processed_rna_bams/{sample}.{genome}.newheader.sam"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = long
	run:
		shell("{params.samtools} view -H {input.rbam} | grep '@HD' > {output}")
		shell("{params.samtools} view -H {input.gbam} | grep '@SQ' >> {output}")
		shell("{params.samtools} view -H {input.rbam} | grep '@RG' >> {output}")
		shell("{params.samtools} view -H {input.rbam} | grep '@PG' >> {output}")

rule rna_bam_reheader:
	input:
		rbam = "processed_rna_bams/{sample}.{genome}.sorted.bam",
		header = "processed_rna_bams/{sample}.{genome}.newheader.sam"
	output:
		"processed_rna_bams/{sample}.{genome}.sorted.reheadered.bam"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.samtools} reheader -P {input.header} {input.rbam} > {output}"

rule index_reheadered_bam:
	input:
		"processed_rna_bams/{sample}.{genome}.sorted.reheadered.bam"
	output:
		"processed_rna_bams/{sample}.{genome}.sorted.reheadered.bam.bai"
	params:
		samtools = samtools_path,
		threads = 4,
		mem = 16,
		t = very_short
	shell:
		"{params.samtools} index {input}"

rule gatk_gvcf_per_scaff_rna:
	input:
		ref = "new_reference/{genome}.fasta",
		bam = "processed_rna_bams/{sample}.{genome}.sorted.reheadered.bam",
		bai = "processed_rna_bams/{sample}.{genome}.sorted.reheadered.bam.bai"
	output:
		"gvcfs_rna/{sample}.{genome}.{scaffold}.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		scaff = "{scaffold}",
		threads = 4,
		mem = 16,
		t = very_long
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""HaplotypeCaller -R {input.ref} -I {input.bam} -L {params.scaff} """
		"""-ERC GVCF --do-not-run-physical-phasing -O {output}"""

rule gatk_combinegvcfs_per_scaffold_rna:
	input:
		ref = "new_reference/{genome}.fasta",
		gvcfs = expand(
			"gvcfs_rna/{sample}.{{genome}}.{{scaffold}}.g.vcf.gz", sample=rna)
	output:
		"combined_gvcfs_rna/{genome}.{scaffold}.gatk.combinegvcf.g.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 4,
		mem = 16,
		t = very_long
	run:
		variant_files = []
		for i in input.gvcfs:
			variant_files.append("--variant " + i)
		variant_files = " ".join(variant_files)
		shell(
			"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
			"""CombineGVCFs -R {input.ref} {variant_files} -O {output}""")

rule gatk_genotypegvcf_per_scaffold_rna:
	input:
		ref = "new_reference/{genome}.fasta",
		gvcf = "combined_gvcfs_rna/{genome}.{scaffold}.gatk.combinegvcf.g.vcf.gz"
	output:
		"genotyped_vcfs_rna/{genome}.{scaffold}.gatk.called.raw.vcf.gz"
	params:
		temp_dir = temp_directory,
		gatk = gatk_path,
		threads = 4,
		mem = 16,
		t = very_long
	shell:
		"""{params.gatk} --java-options "-Xmx15g -Djava.io.tmpdir={params.temp_dir}" """
		"""GenotypeGVCFs -R {input.ref} -V {input.gvcf} -O {output}"""

rule concatenate_split_vcfs_rna:
	input:
		vcf = lambda wildcards: expand(
			"genotyped_vcfs_rna/{gen}.{scaffold}.gatk.called.raw.vcf.gz",
			gen=wildcards.genome,
			scaffold=config["100kb_scaffolds"])
	output:
		"combined_vcfs_rna/combined.{genome}.raw.vcf.gz"
	params:
		bcftools = bcftools_path,
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"{params.bcftools} concat -O z -o {output} {input.vcf}"

rule index_concatenated_vcf_rna:
	input:
		"combined_vcfs_rna/combined.{genome}.raw.vcf.gz"
	output:
		"combined_vcfs_rna/combined.{genome}.raw.vcf.gz.tbi"
	params:
		tabix = tabix_path,
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"{params.tabix} -p vcf {input}"

rule filter_concatenated_vcf_rna:
	input:
		vcf = "combined_vcfs_rna/combined.{genome}.raw.vcf.gz",
		idx = "combined_vcfs_rna/combined.{genome}.raw.vcf.gz.tbi"
	output:
		"combined_vcfs_rna/combined.{genome}.filtered.vcf.gz"
	params:
		bgzip = bgzip_path,
		bcftools = bcftools_path,
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"{params.bcftools} filter -i "
		"'QUAL >= 30 && MQ >= 30 && QD > 2' {input.vcf} | "
		"{params.bcftools} filter -i 'FMT/DP >= 10 & FMT/GQ >= 30' -S . - | "
		"{params.bgzip} > {output}"

rule index_filtered_vcf_rna:
	input:
		"combined_vcfs_rna/combined.{genome}.filtered.vcf.gz"
	output:
		"combined_vcfs_rna/combined.{genome}.filtered.vcf.gz.tbi"
	params:
		tabix = tabix_path,
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"{params.tabix} -p vcf {input}"

rule calc_het_rate_rna:
	input:
		"combined_vcfs_rna/combined.{genome}.filtered.vcf.gz"
	output:
		"results/{genome}.rna.het_rate.txt"
	params:
		suf = "rna",
		threads = 4,
		mem = 16,
		t = long
	shell:
		"python scripts/Calc_het_vcf.py "
		"--vcf {input} "
		"--sexes misc/sample_sexes.txt "
		"--min_sites 5 "
		"--min_ind 1 "
		"--output_file {output} "
		"--suffix {params.suf}"

rule calc_het_rate_dna:
	input:
		"combined_vcfs/combined.{genome}.filtered.vcf.gz"
	output:
		"results/{genome}.dna.het_rate.txt"
	params:
		suf = "dna",
		threads = 4,
		mem = 16,
		t = long
	shell:
		"python scripts/Calc_het_vcf.py "
		"--vcf {input} "
		"--sexes misc/sample_sexes.txt "
		"--min_sites 5 "
		"--min_ind 1 "
		"--output_file {output} "
		"--suffix {params.suf}"

rule compile_stringtie_results:
	input:
		fai = "new_reference/{assembly}.fasta.fai",
		ctabs = lambda wildcards: expand(
			"stringtie_gtfs_{strat}/{sample}_{genome}/t_data.ctab",
			genome=wildcards.assembly,
			strat=wildcards.strategy,
			sample=rna)
	output:
		"results/{assembly}.{strategy}.stringtie_compiled.txt"
	params:
		strat = "{strategy}",
		threads = 4,
		mem = 16,
		t = long
	run:
		ctab_sexes = []
		for i in input.ctabs:
			i_split = i.split("/")[1]
			sample_id = "{}_{}_{}".format(
				i_split.split("_")[0], i_split.split("_")[1], i_split.split("_")[2])
			ctab_sexes.append(config["sexes"][sample_id])
		shell(
			"python scripts/Compile_stringtie_results.py --fai {input.fai} "
			"--output_file {output} --input_files {input.ctabs} "
			"--sex {ctab_sexes} --suffix {params.strat}")

rule compile_stringtie_results_per_transcript:
	input:
		ctabs = lambda wildcards: expand(
			"stringtie_gtfs_{strat}/{sample}_{genome}/t_data.ctab",
			genome=wildcards.assembly,
			strat=wildcards.strategy,
			sample=rna)
	output:
		"results/{assembly}.{strategy}.stringtie_compiled_per_transcript.txt"
	params:
		strat = "{strategy}",
		threads = 4,
		mem = 16,
		t = long
	run:
		ctab_sexes = []
		for i in input.ctabs:
			i_split = i.split("/")[1]
			sample_id = "{}_{}_{}".format(
				i_split.split("_")[0], i_split.split("_")[1], i_split.split("_")[2])
			ctab_sexes.append(config["sexes"][sample_id])
		shell(
			"python scripts/Compile_stringtie_per_transcript.py "
			"--output_file {output} --input_files {input.ctabs} "
			"--sex {ctab_sexes} --suffix {params.strat}")

rule compile_chrom_stats:
	input:
		stats = "xyalign_analyses/{assembly}/results/{assembly}_chrom_stats_count.txt"
	output:
		df = "results/{assembly}.chromstats_compiled.txt",
		html = "results/{assembly}.chromstats_compiled.html",
		html_no_nan = "results/{assembly}.chromstats_compiled.html_nonan.html",
		html_cutoff = "results/{assembly}.chromstats_compiled.html_cutoff.html",
		plot = "results/{assembly}.chromstats_compiled.png",
	params:
		cov_cutoff = "0.95",
		males = lambda wildcards: expand(
			"{sample}.{genome}.mkdup.sorted.bam",
			sample=[x for x in dna if config["sexes"][x] == "male"],
			genome=[wildcards.assembly]),
		females = lambda wildcards: expand(
			"{sample}.{genome}.mkdup.sorted.bam",
			sample=[x for x in dna if config["sexes"][x] == "female"],
			genome=[wildcards.assembly]),
		threads = 4,
		mem = 16,
		t = long
	shell:
		"python scripts/Compile_chromstats_results.py "
		"--input_file {input.stats} "
		"--male_list {params.males} "
		"--female_list {params.females} "
		"--output_dataframe {output.df} "
		"--output_html {output.html} "
		"--output_html_no_nan {output.html_no_nan} "
		"--output_html_cutoff_df {output.html_cutoff} "
		"--coverage_cutoff {params.cov_cutoff} "
		"--plot_title {output.plot}"

rule combine_into_big_dataframe:
	input:
		expr = "results/{assembly}.{strategy}.stringtie_compiled.txt",
		cov = "results/{assembly}.chromstats_compiled.txt",
		het_dna = "results/{assembly}.dna.het_rate.txt",
		het_rna = "results/{assembly}.rna.het_rate.txt"
	output:
		txt = "results/all_compiled.{assembly}.{strategy}.txt"
	params:
		html = "results/all_compiled.{assembly}.{strategy}.html",
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"python scripts/Compile_big_dataframe.py "
		"--input_files {input.cov} {input.expr} {input.het_dna} {input.het_rna} "
		"--output_html_df {params.html} "
		"--output_csv_df {output.txt}"

rule get_komodo:
	output:
		"komodo/komodov1.fa"
	params:
		web_address = lambda wildcards: config["komodo_paths"]["komodov1"],
		initial_output = "komodo/komodov1.fa.gz",
		threads = 1,
		mem = 4,
		t = very_short
	run:
		shell("wget {params.web_address} -O {params.initial_output}")
		shell("gunzip {params.initial_output}")

rule extract_218:
	input:
		ref = "new_reference/{assembly}.fasta",
		fai = "new_reference/{assembly}.fasta.fai"
	output:
		"komodo/{assembly}.scaff218.fasta"
	params:
		scaffold = "218",
		samtools = samtools_path,
		threads = 2,
		mem = 8,
		t = medium
	shell:
		"{params.samtools} faidx {input.ref} {params.scaffold} > {output}"

rule create_lastdb:
	input:
		target = "komodo/komodov1.fa"
	output:
		"komodo/komodoDb.tis"
	params:
		lastdb = lastdb_path,
		database_prefix = "komodo/komodoDb",
		threads = 4,
		mem = 16,
		t = long
	shell:
		"{params.lastdb} {params.database_prefix} {input.target}"

rule run_lastal:
	input:
		t = "komodo/komodoDb.tis",
		q = "komodo/{assembly}.scaff218.fasta"
	output:
		"komodo/komodo_scaff218_{assembly}_align.maf"
	params:
		lastal = lastal_path,
		database_prefix = "komodoDb",
		threads = 8,
		mem = 32,
		t = very_long
	shell:
		"{params.lastal} -a 400 -b 30 -e 4500 {params.database_prefix} {input.q} > {output}"

rule find_par:
	input:
		males = lambda wildcards: expand(
			"xyalign_analyses/{sample}.{genome}/bed/{sample}.{genome}_full_dataframe_depth_mapq_preprocessing.csv",
			genome=wildcards.genome,
			sample=[x for x in dna if config["sexes"][x] == "male"]),
		females = lambda wildcards: expand(
			"xyalign_analyses/{sample}.{genome}/bed/{sample}.{genome}_full_dataframe_depth_mapq_preprocessing.csv",
			genome=wildcards.genome,
			sample=[x for x in dna if config["sexes"][x] == "male"])
	output:
		"par_results/scaffold{scaff}_{genome}.txt"
	params:
		scaffold = "{scaff}",
		threads = 4,
		mem = 16,
		t = medium
	shell:
		"python scripts/Find_par.py --male_list {input.males} "
		"--female_list {input.females} --scaffold {params.scaffold} --output {output}"
