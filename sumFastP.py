#!/home/jiapengc/mambaforge/bin/python3

import pandas as pd
import json

# Read the samplesheet.tsv file into a DataFrame
# Set the 'Sample' column as the index
metaDF = pd.read_csv("samplesheet.tsv", sep="\t", index_col="Sample")

# Initialize an empty DataFrame with the desired columns for metrics
df = pd.DataFrame(columns=[
    "Sample",                # Sample ID
    "bf_total_reads",        # Before filtering: total number of reads
    "bf_total_bases",        # Before filtering: total number of bases
    "bf_q20_bases",          # Before filtering: number of bases with quality score >= 20
    "bf_q30_bases",          # Before filtering: number of bases with quality score >= 30
    "bf_q20_rate",           # Before filtering: rate of bases with quality score >= 20
    "bf_q30_rate",           # Before filtering: rate of bases with quality score >= 30
    "bf_read1_mean_length",  # Before filtering: mean length of Read 1
    "bf_read2_mean_length",  # Before filtering: mean length of Read 2
    "bf_gc_content",         # Before filtering: GC content percentage
    "af_total_reads",        # After filtering: total number of reads
    "af_total_bases",        # After filtering: total number of bases
    "af_q20_bases",          # After filtering: number of bases with quality score >= 20
    "af_q30_bases",          # After filtering: number of bases with quality score >= 30
    "af_q20_rate",           # After filtering: rate of bases with quality score >= 20
    "af_q30_rate",           # After filtering: rate of bases with quality score >= 30
    "af_read1_mean_length",  # After filtering: mean length of Read 1
    "af_read2_mean_length",  # After filtering: mean length of Read 2
    "af_gc_content",         # After filtering: GC content percentage
    "duplication_rate",      # Duplication rate
    "insert_size_peak",      # Insert size peak
    "insert_size_unknown"    # Insert size unknown
])

# Iterate through each sample in the metaDF DataFrame
for s in metaDF.index:
    # Construct the path to the JSON file for the current sample
    dirf = "fastp/" + s + ".json"
    
    # Open and read the JSON file
    with open(dirf) as f:
        jdict = json.load(f)
    
    # Extract the 'before_filtering' and 'after_filtering' metrics from JSON data
    bf_metrics = jdict['summary']['before_filtering']
    af_metrics = jdict['summary']['after_filtering']
    duplication_rate = jdict['duplication'].get('rate', 'N/A')  # Add duplication rate
    insert_size_peak = jdict['insert_size'].get('peak', 'N/A')  # Add insert size peak
    insert_size_unknown = jdict['insert_size'].get('unknown', 'N/A')  # Add insert size unknown
    
    # Create a list of metric values for the current sample, both before and after filtering
    mlist = [
        s,                                   # Sample ID
        str(bf_metrics['total_reads']),      # Before filtering: total reads
        str(bf_metrics['total_bases']),      # Before filtering: total bases
        str(bf_metrics['q20_bases']),        # Before filtering: Q20 bases
        str(bf_metrics['q30_bases']),        # Before filtering: Q30 bases
        str(bf_metrics['q20_rate']),         # Before filtering: Q20 rate
        str(bf_metrics['q30_rate']),         # Before filtering: Q30 rate
        str(bf_metrics['read1_mean_length']),# Before filtering: Read 1 mean length
        str(bf_metrics['read2_mean_length']),# Before filtering: Read 2 mean length
        str(bf_metrics['gc_content']),       # Before filtering: GC content
        str(af_metrics['total_reads']),      # After filtering: total reads
        str(af_metrics['total_bases']),      # After filtering: total bases
        str(af_metrics['q20_bases']),        # After filtering: Q20 bases
        str(af_metrics['q30_bases']),        # After filtering: Q30 bases
        str(af_metrics['q20_rate']),         # After filtering: Q20 rate
        str(af_metrics['q30_rate']),         # After filtering: Q30 rate
        str(af_metrics['read1_mean_length']),# After filtering: Read 1 mean length
        str(af_metrics['read2_mean_length']),# After filtering: Read 2 mean length
        str(af_metrics['gc_content']),       # After filtering: GC content
        str(duplication_rate),               # Duplication rate
        str(insert_size_peak),               # Insert size peak
        str(insert_size_unknown)             # Insert size unknown
    ]
    
    # Convert the list to a DataFrame row and concat it to the main DataFrame
    new_row = pd.DataFrame([mlist], columns=df.columns)
    df = pd.concat([df, new_row], ignore_index=True)

# Save the final DataFrame with all metrics to a TSV file
df.to_csv("Sequencing.metrics.tsv", sep='\t', index=False)
