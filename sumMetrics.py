import pandas as pd
import json

metaDF = pd.read_csv("samplesheet.tsv", sep = "\t", index_col="Sample")

df = pd.DataFrame(columns=["Sample", "raw_reads", "clean_reads", "mapped2HOMD", "mapped_perc"])
for s in metaDF.index:

    # read fastp json
    fpjs = "fastp/" + s + ".json"
    f = open(fpjs)
    fpd = json.load(f)
    raw_reads = fpd['summary']['before_filtering']['total_reads']

    # read map2HOMD json
    hojs = "map2HOMD/" + s + ".json"
    f = open(hojs)
    hopd = json.load(f)
    clean_reads = hopd['general']['reads']['total'] # reads used for eHOMD mapping are trimmed
    mapped2HOMD = hopd['general']['reads']['mapped']['0']

    mapped_perc = mapped2HOMD * 100 / clean_reads


    mlist = [s,
              str(raw_reads),
              str(clean_reads),
              str(mapped2HOMD),
              str(mapped_perc)
              ]

    df = df.append(pd.Series(mlist, index=df.columns), ignore_index=True)




df.to_csv("metrics.tsv", sep='\t', index=False)