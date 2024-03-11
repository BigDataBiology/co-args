import pronto
import matplotlib.pyplot as plt
import polars as pl
import polars.selectors as cs
import numpy as np
import lzma
import seaborn as sns
from scipy import cluster
from collections import defaultdict, Counter

import argnorm.normalizers
from load_data import load_data, ALL_TOOLS

HABITAT = 'human gut'
tool = 'rgi_strict'
tool = 'deeparg'

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

%matplotlib qt

fig, ax = plt.subplots()
Z = cluster.hierarchy.linkage(counts, method='average', optimal_ordering=True)
R = cluster.hierarchy.dendrogram(Z, ax=ax, labels=interesting, leaf_rotation=90, orientation='left')
leaves = R['leaves']
ax.clear()
sns.heatmap(counts[leaves].T[leaves] + 1, ax=ax, cmap='viridis')
ax.set_title(f'Gene co-occurrence in {tool} ({HABITAT})')
fig.tight_layout()
interesting_reordered = [interesting[i] for i in leaves]

aro = pronto.Ontology('aro.obo')
confers_resistance_to_drug_class = aro.get_relationship('confers_resistance_to_drug_class')
confers_resistance_to_antibiotic = aro.get_relationship('confers_resistance_to_antibiotic')

AB_MOLECULE_ARO = 'ARO:1000003'
assert aro[AB_MOLECULE_ARO].name == 'antibiotic molecule'
Ab_classes = set(aro[AB_MOLECULE_ARO].subclasses(1, with_self=False))
name2term = {}
for k in aro.keys():
    if k.startswith('ARO:'):
        name2term[aro[k].name] = k

terms = []
for g in interesting_reordered[:45]:
    if g in name2term:
        terms.append(aro[name2term[g]])

norm = argnorm.normalizers.DeepARGNormalizer()
aro_map_table = norm.get_aro_mapping_table()
aro_map_table['GeneFamily'] = aro_map_table['Original ID'].str.split('|').str[4].str.upper()

cluster = interesting_reordered[:45]

by_family = defaultdict(set)
for _,row in aro_map_table[['GeneFamily', 'ARO']].iterrows():
    by_family[row.GeneFamily].add(row.ARO)

for c in cluster:
    print(c)
    cur_ab_classes = set()
    all_superclasses = set()
    for a in by_family[c]:
        all_superclasses.update(aro[a].superclasses())
    for p in all_superclasses:
        for ab in p.relationships.get(confers_resistance_to_antibiotic, []):
            cur_ab_classes.update( set(ab.superclasses()) & Ab_classes )
        for dc in p.relationships.get(confers_resistance_to_drug_class, []):
            cur_ab_classes.update( set(dc.superclasses()) & Ab_classes )

    for ab in cur_ab_classes:
        print(f'  {ab.name} ({ab.id})')
    print(f'  superclasses:')
    for s in aro[a].superclasses():
        print(f'         {s.name} ({s.id})')

