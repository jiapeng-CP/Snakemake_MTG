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
                expand("map2human/{sample}.bam", sample=SAMPLES),
                #expand("MetaPhlan/{sample}.txt", sample=SAMPLES),
                #expand("map2HOMD/{sample}.json", sample=SAMPLES),



rule fastp:
	input:
		r1 = get_fastq1,
		r2 = get_fastq2

	output:
		r1 = "fastp/{sample}.r1.fq.gz",
		r2 = "fastp/{sample}.r2.fq.gz",
		json = "fastp/{sample}.json",
		html = "fastp/{sample}.html"

	log: "logs/{sample}.fastp.log"
	threads: 8
	shell:
		"mkdir -p fastp \n"
		"/home/jiapengc/.conda/envs/QC/bin/fastp --in1 {input.r1} "
		"--in2 {input.r2} "
		"--out1 {output.r1} "
		"--out2 {output.r2} "
		"--json {output.json} "
		"--html {output.html} "
		"--thread 8"
        

rule map2human: #bowtie2 mapping, sam2bam, bamstat
	input:
		r1 = "fastp/{sample}.r1.fq.gz",
		r2 = "fastp/{sample}.r2.fq.gz"
	output:
		bam = "map2human/{sample}.bam",
		json = "map2human/{sample}.bamstat.json"
	threads: 8
	shell:
		"mkdir -p map2human \n"
		"/home/jiapengc/.conda/envs/biobakery3/bin/bowtie2 -x /data/databases/human/GRCh38_latest_genomic.fna "
		"-1 {input.r1} -2 {input.r2} "
		"--un-conc-gz map2human/{wildcards.sample}_host_removed.fq.gz "
		"--sensitive --threads 4 | "
		"/home/jiapengc/.conda/envs/QC/bin/samtools view -bS -@ 4 > {output.bam} \n"
		"/home/jiapengc/bin/bamstats --cpu 8 --input {output.bam} > {output.json}"

rule map2HOMD:
	input:
		#below fq1&fq2 are the input files, but it won't match the output file names above
		#fq1 = "map2human/{sample}_host_removed.fq.1.gz",
		#fq2 = "map2human/{sample}_host_removed.fq.2.gz"
		#humanstat = "map2human/{sample}.bamstat.json" # this rule doesn't need the human json, put it here as a proxy
		r1 = "fastp/{sample}.r1.fq.gz", # use fastp reads instead of human rm reads
		r2 = "fastp/{sample}.r2.fq.gz"
	output:
		bam = "map2HOMD/{sample}.bam",
		json = "map2HOMD/{sample}.json"
	threads: 8
	shell:
		"mkdir -p map2HOMD \n"
		"/home/jiapengc/.conda/envs/biobakery3/bin/bowtie2 -x /home/jiapengc/db/HOMD/ALL_genomes.fna "
		"-1 {input.r1} -2 {input.r2} "
		"--sensitive --threads 4 | "
		"/home/jiapengc/.conda/envs/QC/bin/samtools view -bS -@ 4 > {output.bam} \n"
		"/home/jiapengc/bin/bamstats --cpu 8 --input {output.bam} > {output.json}"


rule catFq:
	input:
		r1 = "map2human/{sample}_host_removed.fq.1.gz",
		r2 = "map2human/{sample}_host_removed.fq.2.gz"

	output:
		r1r2fq = "catFq/{sample}.r1r2.fq"
	log: "logs/catFq.log"
	shell:
		"mkdir -p catFq \n"
		"zcat {input.r1} {input.r2} > {output} 2> {log}"


rule MetaPhlan:
	input:
		fq = "catFq/{sample}.r1r2.fq"
	output:
		profiletxt = "MetaPhlan/{sample}.txt"
	log: "logs/{sample}.MetaPhlan.log"
	threads: 8
	shell:
		"mkdir -p MetaPhlan \n"
		"metaphlan --nproc 8 --offline --input_type fastq --output_file {output.profiletxt} {input.fq} > {log} 2>&1"


rule humann: # conda activate /home/artemisl/.conda/envs/biobakery
	input:
		fq = "catFq/{sample}.r1r2.fq"
	log: "logs/{sample}.humann.stdouterr.log"
	threads: 8
	output:
		ofolder = "HumanN/{sample}"
	shell:
		"mkdir -p HumanN \n"
		"/home/artemisl/.conda/envs/biobakery/bin/humann "
		"--metaphlan-options=\"--offline\" --threads 8 "
		"--nucleotide-database /home/jiapengc/db/humannDB/chocophlan "
		"--protein-database /home/jiapengc/db/humannDB/uniref "
		"--input {input.fq} --output {output.ofolder} > {log} 2>&1"
