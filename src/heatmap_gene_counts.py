import matplotlib.pyplot as plt
import polars as pl
import polars.selectors as cs
import numpy as np
import lzma
import seaborn as sns
from scipy import cluster, spatial

from load_data import load_data, ALL_TOOLS

HABITAT = 'human gut'

for tool in ALL_TOOLS:
    df = load_data(tool, HABITAT)

    freq = df.select(cs.by_dtype(pl.Int64)) \
            .select(pl.col("*") > 0) \
            .sum()

    nr_samples = df.shape[0]
    min_freq = nr_samples // 5
    max_freq = nr_samples - min_freq

    interesting = [col.name for col in
                        freq.select((pl.col("*") > min_freq) & (pl.col("*") < max_freq))
                    if col.all()]
    nr_mid_freq_genes = len(interesting)

    print(f"{tool}: {df.shape[0]} samples x {df.shape[1]} genes ({nr_mid_freq_genes} genes with 20-80% presence)")
    counts = np.zeros((nr_mid_freq_genes, nr_mid_freq_genes), dtype=np.int64)
    for i in range(nr_mid_freq_genes):
        counts[i,i] = df.select(interesting[i]).select(pl.col("*") > 0).sum()[0,0]
        for j in range(i + 1, nr_mid_freq_genes):
            mask = df.select(
                    [interesting[i], interesting[j]]) \
                    .select(pl.col("*") > 0) \
                    .sum_horizontal()
            both = (mask == 2).sum()
            counts[i,j] = both
            counts[j,i] = both

    fig, ax = plt.subplots()

    Z = cluster.hierarchy.linkage(counts, method='average', optimal_ordering=True)
    R = cluster.hierarchy.dendrogram(Z, ax=ax, labels=interesting, leaf_rotation=90, orientation='left')
    leaves = R['leaves']
    ax.clear()
    sns.heatmap(counts[leaves].T[leaves] + 1, ax=ax, cmap='viridis')
    ax.set_title(f'Gene co-occurrence in {tool} ({HABITAT})')
    fig.tight_layout()
    fig.savefig(f'plots/heatmap_{tool}_{HABITAT}.pdf')

