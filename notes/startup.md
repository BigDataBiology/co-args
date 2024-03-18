# Start up project

## Set up conda/repository

```bash
mkdir co-args
cd co-args
git init 
mkdir notes
conda create -n co-args python=3.11 numpy matplotlib jug
conda activate co-args
conda install polars seaborn pytest pronto requests ipython
```

For running RGI, we need a separate environment, which we will call simply `rgi`:

```bash
conda create -n rgi rgi
```

### Packages

- [Polars](https://pola.rs/)
- [Jug](https://jug.rtfd.io/) and we add [ipython](https://ipython.org/) to enable `jug shell`
- [pytest](https://pytest.org/)
- [pronto](https://pronto.readthedocs.io/)
- [requests](https://requests.readthedocs.io/)

## Retrieve data

### Metagenomes

ARG abundance files are available in the [arg-compare repository](https://github.com/BigDataBiology/arg-compare):

```bash
git clone https://github.com/BigDataBiology/arg-compare
mkdir data
cd data
ln -s ../arg-compare/GMGCv1-unigene-abundance-projection/computed/*.tsv.xz .
cd ..
echo data/ >> .gitignore
```

### GMGCv1 Metatada

GMGCv1 metadata is available for download using the browser at
[https://gmgc.embl.de/download.cgi](https://gmgc.embl.de/download.cgi) but also
using [git-annex](https://git-annex.branchable.com/) at
[GMGC10.data](https://git.embl.de/coelho/GMGC10.data):

```bash
git clone https://git.embl.de/coelho/GMGC10.data.git
cd GMGC10.data/metadata
git-annex get GMGC10.sample.meta.tsv.gz
```

### ARO OBO file

You need to obtain the aro.obo file from https://card.mcmaster.ca/download

Because of a [pronto issue](https://github.com/althonos/pronto/issues/140), we need the OBO format file and not the OWL file.
