import pandas as pd

samplesDF = pd.read_csv("samplesheet.tsv", sep="\t", index_col="Sample")
SAMPLES = []
sampleFq1 = {}
sampleFq2 = {}


for x in samplesDF.index:
	SAMPLES.append(x)
	sampleFq1[x] = samplesDF.loc[x,'fq1']
	sampleFq2[x] = samplesDF.loc[x,'fq2']

def get_sample(wildcards):
	return(wildcards.sample)

def get_fastq1(wildcards):
	return(samplesDF.loc[wildcards.sample,'fq1'])

def get_fastq2(wildcards):
	return(samplesDF.loc[wildcards.sample,'fq2'])

rule all:
	input:
		"allSample.metaphlan.txt",
		expand("HumanN/{sample}", sample=SAMPLES)

rule KneadData:
	input: 
		r1 = get_fastq1,
		r2 = get_fastq2
		
	output:
		p1 = "KneadData/{sample}_paired_1.fastq",
		p2 = "KneadData/{sample}_paired_2.fastq",
		um1 = "KneadData/{sample}_unmatched_1.fastq",
		um2 = "KneadData/{sample}_unmatched_2.fastq"#,
		#rep1 = "KneadData/{sample}.repeats.removed.1.fastq",
		#rep2 = "KneadData/{sample}.repeats.removed.2.fastq",
		#tr1 = "KneadData/{sample}.trimmed.1.fastq",
		#tr2 = "KneadData/{sample}.trimmed.2.fastq",
		
		
	log: "logs/{sample}.KneadData.log"
	threads: 8
	params:
		s = get_sample
	shell:
		"mkdir -p KneadData \n"
		"kneaddata --output-prefix {params.s} -i1 {input.r1} -i2 {input.r2} -o KneadData "
		"--threads 8 "
		"--reference-db /data/databases/human/ " 
		"--trimmomatic=/data/Microbiome/RefData/Metagenomics/biobakery_workflows_databases/Trimmomatic-0.39/ "
		"--trimmomatic-options ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:TRUE "
		"--sequencer-source TruSeq3 "
		"--trf /home/artemisl/.conda/pkgs/trf-4.09.1-hec16e2b_2/bin/ "
		"--log {log} "
		#https://www.reddit.com/r/bioinformatics/comments/ee6x7h/new_line_breaks_for_long_commands_in_snakemake/


rule catFq:
	input:
		p1 = "KneadData/{sample}_paired_1.fastq",
		p2 = "KneadData/{sample}_paired_2.fastq",
		um1 = "KneadData/{sample}_unmatched_1.fastq",
		um2 = "KneadData/{sample}_unmatched_2.fastq"

	output:
		all4fq = "catFq/{sample}.all4.fq"
	log: "logs/catFq.log"
	shell:
		"mkdir -p catFq \n"
		"cat {input.p1} {input.p2} {input.um1} {input.um2} > {output} 2> {log}"


rule MetaPhlan:
	input:
		fq = "catFq/{sample}.all4.fq"
	output:
		profiletxt = "MetaPhlan/{sample}.txt"
	log: "logs/{sample}.MetaPhlan.log"
	threads: 8
	shell:
		"mkdir -p MetaPhlan \n"
		"metaphlan --nproc 8 -t rel_ab_w_read_stats --input_type fastq --output_file {output.profiletxt} {input.fq} > {log} 2>&1"

rule MetaPhlanMerge:
	input:
		expand("MetaPhlan/{sample}.txt", sample=SAMPLES)
	output:
		"allSample.metaphlan.txt"
	shell:
		"merge_metaphlan_tables.py -o allSample.metaphlan.txt MetaPhlan/*.txt"

rule humann: # conda activate /home/artemisl/.conda/envs/biobakery
	input:
		fq = "catFq/{sample}.all4.fq"
	log: "logs/{sample}.humann.stdouterr.log"
	threads: 8
	output:
		ofolder = "HumanN/{sample}"
	shell:
		"mkdir -p HumanN \n"
		"/home/artemisl/.conda/envs/biobakery/bin/humann --threads 8 --input {input.fq} --output {output.ofolder} > {log} 2>&1"



