import matplotlib.pyplot as plt
import polars as pl
import polars.selectors as cs
import numpy as np
import lzma
import seaborn as sns

meta = pl.read_csv("../data/GMGC10.sample.meta.tsv.gz", separator='\t')
meta = meta.filter(pl.col("habitat") == "human gut") \
        .filter(pl.col("insertsHQ") > 2_000_00)

TOOLS = [
    'abricate_ARG_ANNOT',
    'abricate_CARD',
    'abricate_MEGARES',
    'abricate_NCBI',
    'abricate_RESFINDER_FG2',
    'abricate_RESFINDER',
    'deeparg',
    'rgi',
    'rgi_strict',
    ]

fig,axes = plt.subplots(3,3, sharex=True, sharey=True)

for i,tool in enumerate(TOOLS):
    ax = axes[i//3, i%3]
    with lzma.open(f"../data/{tool}.projected.normed10m.tsv.xz", 'rb') as f:
        df = pl.read_csv(f, separator='\t')

    if tool == 'deeparg':
        df = df.rename({'':'sample'})
    df = df.filter(pl.col('sample').is_in(meta.select(['sample_id'])))

    freq = df.select(cs.by_dtype(pl.Int64)) \
            .select(pl.col("*") > 0) \
            .sum(axis=0)
    freq = freq.to_numpy().ravel()
    ax.hist(freq, bins=100)
    ax.set_title(tool)
    ax.set_xlabel("Number of samples")
    ax.set_ylabel("Frequency")

for ax in axes.flat: ax.set_yscale('log')
sns.despine(fig, trim=True)
fig.tight_layout()
