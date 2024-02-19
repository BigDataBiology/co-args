ALL_TOOLS = [
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


def load_data(tool, habitat, min_nr_inserts=2_000_000):
    '''
    Load data for a given ARG annotation tool and habitat

    Parameters
    ----------
    tool : str
        ARG annotation tool to load data for
    habitat : str
        Habitat to filter for
    min_nr_inserts : int
        Minimum number of inserts for a sample to be included

    Returns
    -------
    df : polars.DataFrame
    '''
    import polars as pl
    import lzma
    meta = pl.read_csv("./data/GMGC10.sample.meta.tsv.gz", separator='\t')
    meta = meta.filter(pl.col("habitat") == habitat) \
            .filter(pl.col("insertsHQ") > min_nr_inserts)

    with lzma.open(f"./data/{tool}.projected.normed10m.tsv.xz", 'rb') as f:
        df = pl.read_csv(f, separator='\t')
    if tool == 'deeparg':
        df = df.rename({'':'sample'})
    df = df.filter(pl.col('sample').is_in(meta.select(['sample_id'])))
    return df

