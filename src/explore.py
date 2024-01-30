import matplotlib.pyplot as plt
import polars as pl
import polars.selectors as cs
import numpy as np
import lzma
import seaborn as sns

df = pl.read_csv(lzma.open("../data/deeparg-normed10m.tsv.xz", 'rb'), separator='\t')

freq = df.select(cs.by_dtype(pl.Int64)) \
        .select(pl.col("*") > 0) \
        .sum(axis=0)
freq = freq.to_numpy().ravel()

fig,ax = plt.subplots()
ax.clear()
ax.hist(freq, bins=100)
ax.set_xlabel("Number of samples")
ax.set_ylabel("Frequency")
sns.despine(fig, trim=True)
