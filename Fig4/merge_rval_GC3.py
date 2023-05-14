#script purpose: merge the csv containing the R values with the csv containing the GC3 values

import pandas as pd

columns = ["acc", "rvals", "pvals"]
rvals_df = pd.read_csv("correlation_df.csv", usecols = columns)

GC3_df = pd.read_csv("GC3_content.csv")

GC3_df_dict = GC3_df.set_index("accession").transpose().to_dict(orient = "dict")
rvals_df_dict = rvals_df.set_index("acc").transpose().to_dict(orient = "dict")

df_GC3 = pd.DataFrame(GC3_df_dict)
df_rvals = pd.DataFrame(rvals_df_dict)

df_comb = pd.concat([df_GC3, df_rvals], axis=0)
df_comb = df_comb.transpose()

df_comb.to_csv("rval_vs_GC3.csv")