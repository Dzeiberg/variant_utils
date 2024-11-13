import urllib.request
from pathlib import Path
import pandas as pd

def get_gene_info(geneSymbol : str, cache_dir = ".cache", summary_file:str="MANE.GRCh38.v1.4.summary.txt.gz") -> dict:
    """
    Get information about a gene from the MANE summary file

    Parameters
    ----------
    geneSymbol : str

    Optional Parameters
    -------------------
    mane_summary_file : str
        The path to the MANE summary file or the path to which it should be saved

    Returns
    -------
    dict
        NCBI_GeneID
        Ensembl_Gene
        HGNC_ID
        symbol
        name
        RefSeq_nuc
        RefSeq_prot
        Ensembl_nuc
        Ensembl_prot
        MANE_status
        GRCh38_chr
        chr_start
        chr_end
        chr_strand
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(exist_ok=True)
    mane_summary_filepath = cache_dir / summary_file
    if not mane_summary_filepath.exists():
        urllib.request.urlretrieve(f"https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/{summary_file}",mane_summary_filepath)
    mane_summary = pd.read_csv(mane_summary_filepath, sep="\t",compression='gzip')
    record = mane_summary[mane_summary.symbol == geneSymbol]
    if record.empty:
        raise ValueError(f"Gene {geneSymbol} not found in MANE summary file")
    record = record.iloc[0]
    CHROM = record['GRCh38_chr']
    chrom = str(int(CHROM.split(".")[0].replace("NC_","")))
    if chrom == "23":
        chrom = "X"
    elif chrom == "24":
        chrom = "Y"
    record['CHROM'] = chrom
    return record