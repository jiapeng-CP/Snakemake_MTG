# Metagenomics pipeline with Snakemake
## Overview
This repository hosts a Snakemake-based metagenomics pipeline designed for comprehensive analysis of microbial communities from sequencing data. The pipeline integrates several established bioinformatics tools to perform read trimming, taxonomy profiling, and functional profiling.


## Pipeline components

The current version of the pipeline, considered a minimal viable product (MVP), incorporates the following tools from the Artemis Bioinformatics suite:


*   Reads trimming tool `FastP`
*   Taxonomy profiling `MetaPhLan`
*   Functional profiling `HumanN`

Future iterations will introduce additional features, including taxa normalization and differential abundance analysis, to enhance the pipeline's capabilities.

## Why snakemake

Snakemake is a workflow management system that facilitates the creation and execution of data analysis pipelines. It is Python-based and functions similarly to GNU Make, but with a focus on bioinformatics workflows. Snakemake automatically resolves dependencies between analysis steps, as illustrated by the following rule examples:


```
rule map_reads:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "results/mapped/{sample}.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "bwa mem {input} | samtools view -b - > {output}"


rule sort_alignments:
    input:
        "results/mapped/{sample}.bam"
    output:
        "results/mapped/{sample}.sorted.bam"
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools sort -o {output} {input}"
```

Compared to pipelines build with Bash scripts, Snakemake can executes tasks for individual samples in parallel,
and also perform independent tasks for one sample in paralell. In addition, Snakemake pipeline is more maintainable and scalable when more tools
are added to the pipeline in the future.


## To run the pipeline
### Step 1 build Snakefile
Copy the provided Snakefile to your working directory:

```
git clone https://github.com/jiapeng-CP/Snakemake_MTG && cp Snakemake_MTG/Snakefile yourwd
```
### Step 2 make sample sheet
Within your working directory, create a tab-separated sample sheet named samplesheet.tsv with three columns: Sample, fq1, fq2. A template can be found at https://github.com/jiapeng-CP/Snakemake_MTG/blob/main/samplesheet.tsv.

### Step 3 pre-execute Snakemake
Generate a visual representation of the workflow steps to be executed:

```
/home/jiapengc/mambaforge/envs/snakemake/bin/snakemake --dag | dot -Tsvg > pipeline.svg
```
The resulting SVG graph, as shown below, illustrates the pipeline steps for two example samples (EM9232022LRM01HP1 and EM9232022LRM24hrUT1), starting with KneadData and proceeding through MetaPhLan and HUMAnN in parallel.

![pipeline (1)](https://github.com/jiapeng-CP/Snakemake_MTG/assets/131789717/7f1a7234-8afd-4af0-bd11-25f4184874ff)


### Step 4 execute Snakemake
Activate Biobakery conda environment
```
conda activate /home/jiapengc/bin/biobakery4
```
Run Snakemake that will parse Snakefile in the current directory
```
/home/jiapengc/mambaforge/envs/snakemake/bin/snakemake --core 32 
```
The pipeline will generate output files in four designated folders corresponding to each analysis step:

*   FastP
*   MetaPhlan
*   HumanN

Run the below command to generate a report html that details resources and time used for the tasks
```
/home/jiapengc/mambaforge/envs/snakemake/bin/snakemake --report workflow.report.html
```

