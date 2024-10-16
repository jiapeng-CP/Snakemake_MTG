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
                "Sequencing.metrics.tsv"



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
        shell:
                "/home/jiapengc/mambaforge/bin/python3 sumFastP.py"

rule multiqc:
	input: expand("fastp/{sample}.json", sample=SAMPLES)
	output: "multiqc_report.html"
	shell: "/home/jiapengc/mambaforge/bin/multiqc fastp/*"

rule map2human: #bowtie2 mapping, sam2bam, bamstat
	input:
		r1 = "fastp/{sample}.r1.fq.gz",
		r2 = "fastp/{sample}.r2.fq.gz"
	output:
		bam = "map2human/{sample}.bam",
		json = "map2human/{sample}.bamstat.json",
		unmap = ["map2human/{sample}_host_removed.fq.1.gz", "map2human/{sample}_host_removed.fq.2.gz"]
	threads: 8
	log: "logs/map2human.{sample}.log"
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
	log: "logs/map2HOMD.{sample}.log"
	shell:
		"mkdir -p map2HOMD \n"
		"bowtie2 -x /home/jiapengc/db/HOMD/V10.1/ALL_genomes.fna "
		"-1 {input.r1} -2 {input.r2} "
		"--sensitive --threads 7 2> {log} | "
		"/home/jiapengc/.conda/envs/QC/bin/samtools view -bS -@ 1 > {output.bam} \n"
		"/home/jiapengc/bin/bamstats --cpu 8 --input {output.bam} > {output.json} 2>> {log}"




rule MetaPhlan:
	input:
		r1 = "map2human/{sample}_host_removed.fq.1.gz",
		r2 = "map2human/{sample}_host_removed.fq.2.gz"
	output:
		r1r2fq = "catFq_tmp/{sample}.r1r2.fq",
		profiletxt = "MetaPhlan/{sample}.txt"
	log: "logs/MetaPhlan.{sample}.log"
	threads: 8
	shell:
		"mkdir -p catFq_tmp \n"
		"zcat {input.r1} {input.r2} > {output.r1r2fq} && rm {input.r1} {input.r2} \n"
		"mkdir -p MetaPhlan \n"
		"/home/jiapengc/.conda/envs/biobakery/bin/metaphlan -t rel_ab_w_read_stats --nproc {threads} --offline --input_type fastq --add_viruses --no_map "
		"--bowtie2db /home/jiapengc/db/metaphlan_db --index mpa_vJun23_CHOCOPhlAnSGB_202403 "
		"--output_file {output.profiletxt} {output.r1r2fq} > {log} 2>&1"

rule combine_metaphlan:
    input:
        expand("MetaPhlan/{sample}.txt", sample=SAMPLES)
    output:
        "combined_metaphlan_results.tsv"
    log:
        "logs/combine_metaphlan.log"
    shell:
        "/home/jiapengc/mambaforge/bin/python3 combine_metaphlan_results.py > {log} 2>&1"

rule humann: 
	input:
		fq = "catFq_tmp/{sample}.r1r2.fq", # should set fq as intermediate file. this rule rm (update) the fq at then end, will rerun in next snakemake: reason: Input files updated by another job
		metaphlano = "MetaPhlan/{sample}.txt"
	log: "logs/{sample}.humann.stdouterr.log"
	threads: 8
	output:
		ofolder = directory("HumanN/{sample}")
	shell:
		"mkdir -p HumanN \n"
		"humann "
		"--taxonomic-profile {input.metaphlano} "
		#"--metaphlan-options=\"--offline\ " no need with --taxonomic-profile
		"--threads {threads} --remove-temp-output "
		"--nucleotide-database /home/jiapengc/db/humannDB/chocophlan "
		"--protein-database /home/jiapengc/db/humannDB/uniref "
		"--input {input.fq} --output {output.ofolder} > {log} 2>&1 && "
		"rm {input.fq}"
