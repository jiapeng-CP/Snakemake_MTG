import pandas as pd
import json

metaDF = pd.read_csv("samplesheet.tsv", sep = "\t", index_col="Sample")

df = pd.DataFrame(columns=["Sample", "raw_reads", "clean_reads", "mapped2HOMD", "HOMD_perc"])
for s in metaDF.index:

    # read fastp json
    fpjs = "fastp/" + s + ".json"
    with open(fpjs, 'r') as f:
        fpd = json.load(f)
    raw_reads = fpd['summary']['before_filtering']['total_reads']

    # read map2HOMD json
    hojs = "map2HOMD/" + s + ".json"
    with open(hojs) as f:
        hopd = json.load(f)
    clean_reads = hopd['general']['reads']['total'] # reads used for eHOMD mapping are trimmed
    mapped2HOMD = hopd['general']['reads']['mapped']['0']

    HOMD_perc = mapped2HOMD * 100 / clean_reads

    # read map2human json
    hujs = "map2human/" + s + ".bamstat.json"
    with open(hujs) as f:
        hupd = json.load(f)
    mapped2Human = hupd['general']['reads']['mapped']['0']

    Human_perc = mapped2Human * 100 / clean_reads


    mlist = [s,
              str(raw_reads),
              str(clean_reads),
              str(mapped2HOMD),
              str(HOMD_perc),
              str(mapped2Human),
              str(Human_perc)
              ]

    df = df.append(pd.Series(mlist, index=df.columns), ignore_index=True)




df.to_csv("metrics.tsv", sep='\t', index=False)