import pandas as pd

samplesDF = pd.read_csv("samplesheet.tsv", sep="\t", index_col="Sample")
SAMPLES = []
sampleFq1 = {}
sampleFq2 = {}

for x in samplesDF.index:
	SAMPLES.append(x)
	sampleFq1[x] = samplesDF.loc[x,'fq1']
	sampleFq2[x] = samplesDF.loc[x,'fq2']

def get_fastq1(wildcards):
	return(samplesDF.loc[wildcards.sample,'fq1'])

def get_fastq2(wildcards):
	return(samplesDF.loc[wildcards.sample,'fq2'])

rule all:
	input:
		expand("HumanN/{sample}", sample=SAMPLES),
		expand("MetaPhlan/{sample}.txt", sample=SAMPLES),
		"multiqc_report.html",
		"combined_metaphlan_results.tsv",
		"Sequencing.metrics.tsv",
		expand("Kraken2/{sample}_kraken2_output.txt", sample=SAMPLES),
		expand("bracken_reports/{sample}.breport", sample=SAMPLES),
		"kraken2_otu_table.biom"

rule fastp:
	input:
		r1 = get_fastq1,
		r2 = get_fastq2
	output:
		r1 = "fastp/{sample}.r1.fq.gz",
		r2 = "fastp/{sample}.r2.fq.gz",
		json = "fastp/{sample}.json",
		html = "fastp/{sample}.html"
	log: "logs/fastp.{sample}.log"
	threads: 8
	benchmark: "benchmarks/fastp.{sample}.txt"
	shell:
		"mkdir -p fastp \n"
		"/home/jiapengc/.conda/envs/QC/bin/fastp --in1 {input.r1} "
		"--in2 {input.r2} "
		"--out1 {output.r1} "
		"--out2 {output.r2} "
		"--json {output.json} "
		"--html {output.html} "
		"--thread {threads} "
		"2> {log}"

rule FqMetrics:
	input:
		json = expand("fastp/{sample}.json", sample=SAMPLES)
	output:
		"Sequencing.metrics.tsv"
	benchmark: "benchmarks/FqMetrics.txt"
	shell:
		"/home/jiapengc/mambaforge/bin/python3 sumFastP.py"

rule multiqc:
	input: expand("fastp/{sample}.json", sample=SAMPLES)
	output: "multiqc_report.html"
	benchmark: "benchmarks/multiqc.txt"
	shell: "/home/jiapengc/mambaforge/bin/multiqc fastp/*"

rule map2human: #bowtie2 mapping, sam2bam, bamstat
	input:
		r1 = "fastp/{sample}.r1.fq.gz",
		r2 = "fastp/{sample}.r2.fq.gz"
	output:
		bam = "map2human/{sample}.bam",
		json = "map2human/{sample}.bamstat.json",
		temp_unmap = temp(["map2human/{sample}_host_removed.fq.1.gz", "map2human/{sample}_host_removed.fq.2.gz"])
	threads: 8
	log: "logs/map2human.{sample}.log"
	benchmark: "benchmarks/map2human.{sample}.txt"
	shell:
		"mkdir -p map2human \n"
		"bowtie2 -x /data/databases/human/GRCh38_latest_genomic.fna "
		"-1 {input.r1} -2 {input.r2} "
		"--un-conc-gz map2human/{wildcards.sample}_host_removed.fq.gz "
		#"--met-stderr --quiet " # met is not a summary of the alignment e.g. mapping rate etc.
		"--sensitive --threads 4 2> {log} | "
		"/home/jiapengc/.conda/envs/QC/bin/samtools view -bS -@ 4 > {output.bam} \n"
		"/home/jiapengc/bin/bamstats --cpu 8 --input {output.bam} > {output.json} 2>> {log}"

rule map2HOMD:
	input:
		r1 = "fastp/{sample}.r1.fq.gz", # use fastp reads instead of human rm reads
		r2 = "fastp/{sample}.r2.fq.gz"
	output:
		bam = "map2HOMD/{sample}.bam",
		json = "map2HOMD/{sample}.json"
	threads: 8
	log: "logs/map2HOMD.{sample}.log"
	benchmark: "benchmarks/map2HOMD.{sample}.txt"
	shell:
		"mkdir -p map2HOMD \n"
		"bowtie2 -x /home/jiapengc/db/HOMD/V10.1/ALL_genomes.fna "
		"-1 {input.r1} -2 {input.r2} "
		"--sensitive --threads 7 2> {log} | "
		"/home/jiapengc/.conda/envs/QC/bin/samtools view -bS -@ 1 > {output.bam} \n"
		"/home/jiapengc/bin/bamstats --cpu 8 --input {output.bam} > {output.json} 2>> {log}"

rule concat_r1r2:
	input:
		r1 = "map2human/{sample}_host_removed.fq.1.gz",
		r2 = "map2human/{sample}_host_removed.fq.2.gz"
	output:
		r1r2fq = temp("catFq_tmp/{sample}.r1r2.fq")
	benchmark: "benchmarks/concat_r1r2.{sample}.txt"
	shell:
		"""
		mkdir -p catFq_tmp
		zcat {input.r1} {input.r2} > {output.r1r2fq}
		"""

rule MetaPhlan:
	input:
		r1r2fq = "catFq_tmp/{sample}.r1r2.fq"
	output:
		profiletxt = "MetaPhlan/{sample}.txt"
	log: "logs/MetaPhlan.{sample}.log"
	threads: 8
	benchmark: "benchmarks/MetaPhlan.{sample}.txt"
	shell:
		"""
		mkdir -p MetaPhlan
		/home/jiapengc/.conda/envs/biobakery/bin/metaphlan -t rel_ab_w_read_stats --nproc {threads} --offline --input_type fastq --add_viruses --no_map \
		--bowtie2db /home/jiapengc/db/metaphlan_db --index mpa_vJun23_CHOCOPhlAnSGB_202403 \
		--output_file {output.profiletxt} {input.r1r2fq} > {log} 2>&1
		"""

rule combine_metaphlan:
	input:
		expand("MetaPhlan/{sample}.txt", sample=SAMPLES)
	output:
		"combined_metaphlan_results.tsv"
	log:
		"logs/combine_metaphlan.log"
	benchmark: "benchmarks/combine_metaphlan.txt"
	shell:
		"/home/jiapengc/mambaforge/bin/python3 combine_metaphlan_results.py > {log} 2>&1"

rule humann:
	input:
		fq = "catFq_tmp/{sample}.r1r2.fq", 
		metaphlano = "MetaPhlan/{sample}.txt"
	log: "logs/{sample}.humann.stdouterr.log"
	threads: 8
	output:
		ofolder = directory("HumanN/{sample}")
	benchmark: "benchmarks/humann.{sample}.txt"
	shell:
		"""
		mkdir -p HumanN
		humann \
		    --taxonomic-profile {input.metaphlano} \
		    --threads {threads} --remove-temp-output \
		    --nucleotide-database /home/jiapengc/db/humannDB/chocophlan \
		    --protein-database /home/jiapengc/db/humannDB/uniref \
		    --input {input.fq} --output {output.ofolder} > {log} 2>&1
		"""

# kraken2 needs 74.4g mem for 1 sample (db k2_pluspf_20231009)
rule kraken2:
	input:
		r1 = "map2human/{sample}_host_removed.fq.1.gz",
		r2 = "map2human/{sample}_host_removed.fq.2.gz"
	output:
		report = "Kraken2/{sample}_kraken2_report.txt",
		output = "Kraken2/{sample}_kraken2_output.txt"
	log:
		"logs/kraken2_{sample}.log"
	threads: 8
	benchmark: "benchmarks/kraken2.{sample}.txt"
	shell:
		"""
		mkdir -p Kraken2
		/home/jiapengc/miniforge3/envs/kraken2/bin/kraken2 --quick \
			--db /home/jiapengc/db/Kraken2/k2_pluspf_20231009 \
			--threads {threads} \
			--report {output.report} \
			--output {output.output} \
			--paired {input.r1} {input.r2} > {log} 2>&1
		"""

rule bracken:
	input:
		report = "Kraken2/{sample}_kraken2_report.txt"
	output:
		bracken_output = "bracken_outputs/{sample}.bracken",
		bracken_report = "bracken_reports/{sample}.breport"
	log:
		"logs/bracken_{sample}.log"
	params:
		bracken_db = "/home/jiapengc/db/Kraken2/k2_pluspf_20231009",
		read_length = 150,
		level = "S",
		threshold = 10
	shell:
		"""
		mkdir -p bracken_outputs bracken_reports
		/home/jiapengc/miniforge3/envs/kraken2/bin/bracken \
			-d {params.bracken_db} \
			-i {input.report} \
			-r {params.read_length} \
			-l {params.level} \
			-t {params.threshold} \
			-o {output.bracken_output} \
			-w {output.bracken_report} > {log} 2>&1
		"""



rule kraken2biom:
	input:
		expand("Kraken2/{sample}_kraken2_report.txt", sample=SAMPLES)
	output:
		"kraken2_otu_table.biom"
	log:
		"logs/kraken2biom.log"
	shell:
		"""
		mkdir -p logs
		/home/jiapengc/miniforge3/envs/kraken2/bin/kraken-biom {input} --fmt json -o {output} > {log} 2>&1
		"""
