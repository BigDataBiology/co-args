# Stream notes

## 2024-04-19 

We wrote code to download [proGenomes](https://progenomes.embl.de/), organised by specI cluster (a specI cluster is going to be our definition of _species_).

We created a separate conda environment to run rgi (called `rgi`) and use `conda run` to call `rgi`. We then wrote code to run `rgi` on every genome.

## 2024-02-20 (Stream 4)

We generated heatmaps for all the tools.
RGI & deepARG have the most genes in the 20-80% range.

For next time:

- [x] use the antibiotic resistance ontology (ARO) to look at some clusters
- [x] use argNorm to normalize deepARG results to ARO

Eventually:

- Analyse isolate genomes

## 2024-02-13 (Stream 3)

For next time:

- [x] start internal "library" (`load_data.py`)
