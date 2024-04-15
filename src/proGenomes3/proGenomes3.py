from jug import TaskGenerator, bvalue
import os
import atomicwrites

MAX_NR_GENOMES = 10

@TaskGenerator
def download_clustering_table():
    import requests
    os.makedirs('data', exist_ok=True)
    URL = 'https://progenomes.embl.de/data/proGenomes3_specI_clustering.tab.bz2'
    r = requests.get(URL, allow_redirects=True, stream=True)
    oname = 'data/proGenomes3_specI_clustering.tab.bz2'
    with open(oname, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)
    return oname

@TaskGenerator
def parse_clustering_table(oname):
    import bz2
    specIs = {}
    with bz2.open(oname, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue
            specI,genomes = line.strip().split('\t')
            specIs[specI] = genomes.split(';')
    return specIs

def sample_genomes(genomes, max_nr_genomes):
    if len(genomes) <= max_nr_genomes:
        return genomes
    import random
    genomes = genomes[:]
    genomes.sort()
    random.seed(genomes[max_nr_genomes])
    return random.sample(genomes, max_nr_genomes)

@TaskGenerator
def download_genomes(specI, genomes, max_nr_genomes):
    import requests
    os.makedirs(f'data/{specI}', exist_ok=True)
    downloaded = []
    for g in sample_genomes(genomes, max_nr_genomes):
        tax = g.split('.')[0]
        oname = f'data/{specI}/{g}.fna.gz'
        if os.path.exists(oname):
            continue
        url = f'https://progenomes.embl.de/dumpSequence.cgi?p={g}&t=c&a={tax}'
        r = requests.get(url, allow_redirects=True, stream=True)
        with atomicwrites.atomic_write(oname, mode='wb', overwrite=True) as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        downloaded.append(oname)
    return downloaded

@TaskGenerator
def run_rgi_for_all_genomes(genomes):
    import subprocess
    results = []
    for g in genomes:
        rgi_out = g.replace('.fna.gz', '.rgi.out')
        results.append(rgi_out)
        if os.path.exists(rgi_out+'.json'):
            continue
        subprocess.check_call(
            ['conda', 'run', '-n', 'rgi',
             'rgi', 'main',
             '-t', 'contig',
             '-a', 'DIAMOND',
             '-n', '8',
             '--include_loose',
             '--clean',
             '--split_prodigal_jobs',
             '--input_sequence', g,
             '--output_file', rgi_out+'-tmp'])
        os.rename(rgi_out+'-tmp.txt', rgi_out+'.txt')
        # We do not keep the JSON versions because they take up too much diskspace
    return results

@TaskGenerator
def concat_partials(fs: list[str]) -> str:
    import polars as pl
    import polars.selectors as cs
    import pathlib
    import lzma

    partials = []
    for f in fs:
        path = pathlib.Path(f + '.txt')
        genome = path.stem.replace('.rgi.out', '')
        table = pl.read_csv(path, separator='\t', dtypes={'Nudged':bool})
        table = table.with_columns(pl.lit(genome).alias('genome'))
        partials.append(table)
    table = pl.concat(partials)
    oname = path.parent / 'rgi-concat.tsv.xz'
    with lzma.open(oname, 'wt') as f:
        table.write_csv(f, separator='\t')
    return oname


clustering_table = parse_clustering_table(download_clustering_table())
clusters = bvalue(clustering_table)

rgi_results = []
for k,genomes in clusters.items():
    spI = download_genomes(k, genomes, MAX_NR_GENOMES)
    rgi_results.append(concat_partials(run_rgi_for_all_genomes(spI)))


