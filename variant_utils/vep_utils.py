from pathlib import Path
import subprocess
from datetime import datetime
import pandas as pd
from io import StringIO
from pysam import VariantFile

def run_vep(df,assembly,**kwargs):
    write_dir = Path(kwargs.get("write_dir","/tmp"))
    write_dir.mkdir(exist_ok=True)
    vep_input_file = write_dir / f"vep_input_{str(datetime.now()).replace(' ','_').replace('.','_').replace(':','_')}.vcf"
    vep_output_file = write_dir / f"vep_output_{str(datetime.now()).replace(' ','_').replace('.','_').replace(':','_')}.vcf"
    df.loc[:,['CHROM','POS','ID','REF','ALT']].to_csv(vep_input_file,sep='\t',index=False)
    cmd = f"docker run -v /data/dbs/vep_data/:/data -v {write_dir}:/inputs ensemblorg/ensembl-vep vep --input_file /inputs/{vep_input_file.name} --output_file STDOUT --format vcf --vcf --offline --assembly {assembly} --everything --fork 4 --cache --dir_cache /data --fasta /data/fasta/Homo_sapiens.GRCh38.cdna.all.fa.bgz --force_overwrite --no_stats"
    print(f"running '{cmd}'")
    result = subprocess.run(cmd.split(" "), stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    with open(vep_output_file,'w') as f:
        f.write(result.stdout)
    vcf = VariantFile(vep_output_file)
    vep_cols = vcf.header.info['CSQ'].description.split("Format: ")[1].split("|")
    
    vep_results = [line for line in result.stdout.split("\n") if line[:2] != "##"]
    vep_df = pd.read_csv(StringIO('\n'.join(vep_results)),delimiter='\t').dropna(subset='INFO')

    vep_series = vep_df.INFO.apply(lambda r: list(map(lambda s: dict(zip(vep_cols,s.split('|'))),r.split(","))))
    info_df = pd.DataFrame(vep_series,index=vep_df.index).explode('INFO')

    info_df = pd.DataFrame.from_records(info_df.INFO.values,index=info_df.index)
    vep_df = pd.merge(vep_df,info_df,left_index=True,right_index=True,validate='one_to_many')
    vep_df = vep_df.assign(CHROM=vep_df.loc[:,'#CHROM'].astype(str),
                            POS=vep_df.loc[:,'POS'].astype(str),
                            REF=vep_df.loc[:,'REF'].astype(str),
                            ALT=vep_df.loc[:,'ALT'].astype(str))
    return vep_df