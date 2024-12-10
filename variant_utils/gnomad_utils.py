from variant_utils.utils import read_external_config
from pathlib import Path
import subprocess
from datetime import datetime
import pandas as pd
from pysam import VariantFile
from typing import List


def queryGnomAD(assembly, CHROM,START,STOP,HGNC_ID,external_config_filepath,**kwargs):
    """
    Query gnomAD for missense variants in a gene; if assembly is 'GRCh37' gnomAD v2.1.1 is used, otherwise gnomAD v4.1 is used
    
    Dependencies:
    - GATK
    - picard
    - java
    - gnomAD exomes and genomes VCF files

    Parameters:
    -----------
    assembly : str
        The genome assembly to use for the query, either 'GRCh37' or 'GRCh38'
    CHROM : str
        The chromosome for which to query gnomAD
    START : int
        The minimum position in the chromosome for which to query gnomAD
    STOP : int
        The maximum position in the chromosome for which to query gnomAD
    Steps:
    1) Get the chromosomal coordinates of the gene from the MANE GTF file
    2) Use GATK SelectVariants to extract variants in the gene from the gnomAD exomes and genomes VCF files
        2A) Filter for SNPs
        2B) Exclude filtered variants
    3) Use GATK MergeVcfs to combine the exomes and genomes VCF files
    4) Use GATK VariantsToTable to convert the combined VCF file to a TSV file
    5) Manually parse the VEP annotations in the TSV file
    6) Filter for missense variants

    Required args:
    - assembly: str : The genome assembly to use for the query, either 'GRCh37' or 'GRCh38'
    - CHROM: str : The chromosome for which to query gnomAD
    - START: int : The minimum position in the chromosome for which to query gnomAD
    - STOP: int : The maximum position in the chromosome for which to query gnomAD

    Optional kwargs:
    - write_dir: str : Path to the directory where the output files will be written : default "/tmp"

    Returns:
    - missense_df: pd.DataFrame : A DataFrame containing parsed VEP annotations for matched missense variants in gnomAD exomes and genomes
    """
    external_tools = read_external_config(external_config_filepath)
    write_dir = Path(kwargs.get("write_dir","/tmp"))
    write_dir.mkdir(exist_ok=True)
    # java = Path(external_tools.get("java"))
    # picard_filepath = Path(external_tools.get("picard_filepath"))
    # assert picard_filepath.exists(), "picard_filepath does not exist"

    release_version = "v4.1" if assembly == "GRCh38" else "r2.1.1"
    if release_version == "r2.1.1":
        chr = ""
        gnomad_vcf_root = Path(external_tools['gnomad_v2_vcf_root'])
    else:
        chr = "chr"
        gnomad_vcf_root = Path(external_tools['gnomad_v4_vcf_root'])
    assert gnomad_vcf_root.exists(), "gnomad_vcf_root does not exist: {}".format(gnomad_vcf_root)
    gnomAD_exomes_filepath = gnomad_vcf_root / f"exomes/gnomad.exomes.{release_version}.sites.{chr}{CHROM}.vcf.bgz"
    gnomAD_genomes_filepath = gnomad_vcf_root / f"genomes/gnomad.genomes.{release_version}.sites.{chr}{CHROM}.vcf.bgz"
    exomes_output_File = write_dir / f"selectvariants_{str(datetime.now()).replace(' ','_')}.exomes.vcf"
    genomes_output_File = write_dir / f"selectvariants_{str(datetime.now()).replace(' ','_')}.genomes.vcf"
    # cmd = f"{external_tools.get('gatk')} SelectVariants -V {gnomAD_exomes_filepath} -L {chr}{CHROM}:{START}-{STOP} --select-type-to-include SNP --exclude-filtered --output {exomes_output_File}"
    cmd = f"docker run -v {gnomAD_exomes_filepath}:/mnt/exomes.vcf -v {exomes_output_File}:/mnt/exomes_output.vcf broadinstitute/gatk gatk SelectVariants -V /mnt/exomes.vcf -L {chr}{CHROM}:{START}-{STOP} --select-type-to-include SNP --exclude-filtered --output /mnt/exomes_output.vcf"
    subprocess.run(cmd.split(" "))
    # cmd = f"{external_tools.get('gatk')} SelectVariants -V {gnomAD_genomes_filepath} -L {chr}{CHROM}:{START}-{STOP} --select-type-to-include SNP --exclude-filtered --output {genomes_output_File}"
    cmd = f"docker run -v {gnomAD_genomes_filepath}:/mnt/genomes.vcf -v {genomes_output_File}:/mnt/genomes_output.vcf broadinstitute/gatk gatk SelectVariants -V /mnt/genomes.vcf -L {chr}{CHROM}:{START}-{STOP} --select-type-to-include SNP --exclude-filtered --output /mnt/genomes_output.vcf"
    subprocess.run(cmd.split(" "))
    output_File = write_dir / f"combinevariants_{str(datetime.now()).replace(' ','_')}.vcf"
    # cmd = f'{java} -jar {picard_filepath} MergeVcfs I={exomes_output_File} I={genomes_output_File} O={output_File}'
    cmd = f"docker run -v {exomes_output_File}:/mnt/exomes_output.vcf -v {genomes_output_File}:/mnt/genomes_output.vcf -v {output_File}:/mnt/output.vcf broadinstitute/gatk gatk MergeVcfs I=/mnt/exomes_output.vcf I=/mnt/genomes_output.vcf O=/mnt/output.vcf"
    subprocess.run(cmd.split(" "))
    tsvout = str(output_File).replace('.vcf','.tsv')
    # variants2table = f"{external_tools.get('gatk')} VariantsToTable -V {output_File} -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -ASF AC -ASF AF -ASF vep -O {tsvout}"
    variants2table = f"docker run -v {output_File}:/mnt/output.vcf -v {tsvout}:/mnt/output.tsv broadinstitute/gatk gatk VariantsToTable -V /mnt/output.vcf -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -ASF AC -ASF AF -ASF vep -O /mnt/output.tsv"
    subprocess.run(variants2table.split(" "))
    gnomAD_df = pd.read_csv(tsvout,delimiter='\t')
    vep_columns = get_vep_columns_from_vcf_header(output_File)
    vep_df = parse_vep(gnomAD_df,columns=vep_columns)
    gnomAD_df = pd.merge(gnomAD_df,vep_df,left_index=True,right_on='index',validate='one_to_many')
    gene_df = gnomAD_df[gnomAD_df.HGNC_ID == HGNC_ID]
    gene_df = gene_df.assign(CHROM=gene_df.CHROM.astype(str).str.replace("chr",""),
                                POS=gene_df.POS.astype(str),
                                REF=gene_df.REF.astype(str),
                                ALT=gene_df.ALT.astype(str))
    return gene_df

def get_vep_columns_from_vcf_header(vcf_file:str)->list:
    """
    Read a gnomAD vcf file and extract the VEP columns from the header

    Parameters:
    -----------
    vcf_file : str
        The path to the gnomAD VCF file

    Returns:
    --------
    list : A list of VEP columns
    """
    vcf = VariantFile(vcf_file)
    return vcf.header.info['vep'].description.split("Format: ")[1].split("|")
    
def parse_vep(df:pd.DataFrame,columns:List[str])->pd.DataFrame:
    """
    parse the 'vep' column of the gnomAD dataframe into its own dataframe

    Parameters:
    -----------
    df : pd.DataFrame
        The gnomAD dataframe
    columns : list
        The VEP columns
    
    Returns:
    --------
    pd.DataFrame : A DataFrame containing the parsed VEP annotations
    """
    vep_series = df.vep.apply(lambda r: list(map(lambda s: dict(zip(columns,s.split('|'))),r.split(","))))
    vep_df = pd.DataFrame(vep_series,index=df.index).explode('vep')
    vep_df = pd.DataFrame.from_records(vep_df.vep.values,index=vep_df.index).reset_index()
    return vep_df
