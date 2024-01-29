# Start up project

## Set up conda/repository

```bash
mkdir co-args
cd co-args
git init 
mkdir notes
conda create -n co-args python=3.11 numpy matplotlib jug 
conda activate co-args
conda install polars seaborn
```

### Packages

- [Polars](https://pola.rs/)
- [Jug](https://jug.rtfd.io/)


## Retrieve data

### Metagenomes

```bash
git clone https://github.com/BigDataBiology/arg-compare
mkdir data
cd data
ln -s ../arg-compare/GMGCv1-unigene-abundance-projection/computed/deeparg-normed10m.tsv.xz
cd ..
echo data/ >> .gitignore
```


