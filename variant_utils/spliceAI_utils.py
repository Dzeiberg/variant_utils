from variant_utils.utils import read_external_config
from datetime import datetime
from pathlib import Path
import subprocess
import pandas as pd

def querySpliceAI(assembly:str, chrom:str, position_min:int, position_max:int,spliceAIRoot:str|Path,**kwargs):
    """
    Query SpliceAI for spliceAI scores in a region of the genome

    Required args:
    - assembly: str : The genome assembly to use for the query, either 'GRCh37' or 'GRCh38'
    - chrom: str : The chromosome for which to query SpliceAI
    - position_min: int : The minimum position in the chromosome for which to query SpliceAI
    - position_max: int : The maximum position in the chromosome for which to query SpliceAI
    - spliceAIRoot : str|Path : Path to the spliceAI VCF files

    Optional kwargs:
    - cache_dir: str : Path to the directory where the output files will be written : default "/tmp"
    
    """
    spliceAIRoot = Path(spliceAIRoot)
    if assembly == 'GRCh38':
        spliceAI_filepath = spliceAIRoot / "spliceai_scores.raw.snv.hg38.vcf.gz"
    elif assembly == 'GRCh37':
        spliceAI_filepath = spliceAIRoot / "spliceai_scores.raw.snv.hg19.vcf.gz"
    assert spliceAI_filepath.exists(), f"spliceAI_filepath does not exist, {spliceAI_filepath}"
    cache_dir = Path(kwargs.get('cache_dir',"/tmp/"))
    cache_dir.mkdir(exist_ok=True)
    output_filepath = cache_dir / f'splice_ai_query_result.{str(datetime.now()).replace(" ","_")}.vcf'
    cmd = f"docker run -v {spliceAI_filepath.parent}:/mnt -v {output_filepath.parent}:/out broadinstitute/gatk gatk SelectVariants -V /mnt/{spliceAI_filepath.name} -L {chrom}:{max(position_min,1)}-{position_max} --output /out/{output_filepath.name}"
    subprocess.run(cmd.split(" "))
    result_df = pd.read_csv(output_filepath,comment='#',delimiter='\t',header=None,
                    names='CHROM POS ID REF ALT QUAL FILTER INFO'.split(" "),
                        dtype={k : str for k in 'CHROM POS REF ALT'.split(" ")})
    result_df = result_df.assign(spliceAI_score=result_df.INFO.apply(lambda s: max(list(map(float,
                                                                                        s.split("|")[2:6])))))
    if kwargs.get('gene_name',None) is not None:
        result_df = result_df.assign(gene_name=result_df.INFO.str.split("|").str[1])
        result_df = result_df[result_df.gene_name == kwargs['gene_name']]
    return result_df