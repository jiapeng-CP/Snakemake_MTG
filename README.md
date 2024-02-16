# Metagenomics pipeline with Snakemake
## Components of the pipeline

In the current minimal viable product (MVP), I used components from Artemis Bioberkery pipeline including

*   Reads trimming tool `KneadData`
*   Taxonomy profiling `MetaPhLan`
*   Functional profiling `HumanN`

In the near future for the next iterations, more components such as taxa normalization and differential abundance analysis will be added.

## What is snakemake and why snakemake

I use [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) to build the pipeline.
It is a popular workflow management system and a Python-based tool used for creating and executing data analysis workflows.
Like `GNU Make` (used for compilation of software), Snakemake can figure out the dependencies of Bioinformatics steps based on rules.
For example, a rule file like below execute `bwa mem` and then `samtools sort` based on the input and output dependency.

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
Copy my Snakefile to your workding directory

```
cp /home/jiapengc/SNpipeline/Snakefile yourwd
```
### Step 2 make sample sheet
In your workding directory, make a sample sheet with file name `samplesheet.tsv`. It should be tab-seperated file with 3 columns: Sample, fq1, fq2. A demo is `/home/jiapengc/SNpipeline/samplesheet.tsv`

### Step 3 pre-execute Snakemake
In the directory, running
```
/home/jiapengc/mambaforge/envs/snakemake/bin/snakemake --dag | dot -Tsvg > pipeline.svg
```
will produce a svg graph to show the specific steps that will be executed. Below is an example svg you can find at /home/jiapengc/SNpipeline/pipeline.svg.


![pipeline (1)](https://github.com/jiapeng-CP/Snakemake_MTG/assets/131789717/7f1a7234-8afd-4af0-bd11-25f4184874ff)


In the figure, two samples `EM9232022LRM01HP1` and `EM9232022LRM24hrUT1` are analyzed by the Snakemake pipeline, starting with KneadData, concatenate FASTQ files produced by KneadData, then MetaPhlan and HumanN in parallel.
### Step 4 execute Snakemake
Activate Artemis Biobakery conda environment
```
conda activate /home/artemisl/.conda/envs/biobakery
```
Run Snakemake that will parse Snakefile in the current directory
```
/home/jiapengc/mambaforge/envs/snakemake/bin/snakemake --core 32
```
This command will execute the pipeline and it will create 4 folders where it deposites output files of the 4 steps.
*   KneadData
*   catFq
*   MetaPhlan
*   HumanN


## Demo
Test of the MVP pipeline was finished in Hill's Server under /home/jiapengc/SNpipeline.
All necessary results were produced in the sub-folders.


## New components to add

*   Reads normalization for taxa and funtion term comparison across different samples
*   Differential abundance analysis to compare Control vs. Treatment sample groups

