from jug import TaskGenerator, bvalue
import os

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

@TaskGenerator
def download_genomes(specI, genomes):
    import requests
    os.makedirs(f'data/{specI}', exist_ok=True)
    for g in genomes:
        tax = g.split('.')[0]
        oname = f'data/{specI}/{g}.fna.gz'
        if os.path.exists(oname):
            continue
        url = f'https://progenomes.embl.de/dumpSequence.cgi?p={g}&t=c&a={tax}'
        r = requests.get(url, allow_redirects=True, stream=True)
        with open(oname, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return f'data/{specI}'

clustering_table = parse_clustering_table(download_clustering_table())
clusters = bvalue(clustering_table)

for k,genomes in clusters.items():
    download_genomes(k, genomes)

