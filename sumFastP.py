import pandas as pd
import json

metaDF = pd.read_csv("samplesheet.tsv", sep = "\t", index_col="Sample")

df = pd.DataFrame(columns=["Sample", "total_reads"])
for s in metaDF.index:
    dirf = "fastp/" + s + ".json"
    f = open(dirf)
    jdict = json.load(f)
    mlist = [s, str(jdict['summary']['before_filtering']['total_reads'])]
    result = "\t".join(mlist)
    df = df.append(pd.Series(mlist, index=df.columns), ignore_index=True)


df.to_csv("fastp.allsample.tsv", sep='\t', index=False)