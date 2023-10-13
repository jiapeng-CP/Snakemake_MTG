import pandas as pd
import glob

metaDF = pd.read_csv("mtt.metadata.tsv", sep = "\t", index_col="SampleID")
# from run samplesheet, use sed -n '/Sample_ID/,$p' /data/Biosolutions/NextSeqData_2023/231005_NS500274_0085_AH2NNWBGXV/SampleSheet.csv
SAMPLES = []



fq_path = "/home/kristenm/Colgate_NextSeq_2022/221030_NS500274_0073_AH55YKBGXK/fastq_files/"


# Initialize an empty DataFrame
df = pd.DataFrame(columns=["Sample", "fq1", "fq2"])

for x in metaDF.index:
        SAMPLES.append(x)
        r1_pattern = fq_path + x + "*" + "_R1_" + "*"
        r2_pattern = fq_path + x + "*" + "_R2_" + "*"


        r1_path = glob.glob(r1_pattern)
        r2_path = glob.glob(r2_pattern)


        mlist = [x, r1_path[0], r2_path[0]]

        result = "\t".join(mlist)
        df = df.append(pd.Series(mlist, index=df.columns), ignore_index=True)

df.to_csv("testsamplesht.tsv", sep='\t', index=False)
