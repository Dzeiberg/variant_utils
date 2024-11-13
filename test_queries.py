from variant_utils.get_gene_info import get_gene_info
from variant_utils.clinvar_utils import queryClinVarVCF
from variant_utils.gnomad_utils import queryGnomAD
from variant_utils.spliceAI_utils import querySpliceAI
import urllib.request
from pathlib import Path

def test_get_gene_info():
    get_gene_info('BRCA1').to_json(".cache/BRCA1.json")

def test_queryClinVarVCF():
    brca1_info = get_gene_info("BRCA1")
    # set destination to save/reload ClinVar
    cache_dir = Path(".cache")
    cache_dir.mkdir(exist_ok=True)
    # download ClinVar release (e.g., 2018-12-17) if file is not present
    clinvar_filepath = cache_dir / "clinvar_20181217.vcf.gz"
    idx_filepath = cache_dir / "clinvar_20181217.vcf.gz.tbi"
    if not clinvar_filepath.exists():
        urllib.request.urlretrieve(f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20181217.vcf.gz",str(clinvar_filepath))
    if not idx_filepath.exists():
        urllib.request.urlretrieve(f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_20181217.vcf.gz.tbi",str(idx_filepath))

    brca1_clinvar_variants = queryClinVarVCF(str(clinvar_filepath), brca1_info.CHROM, brca1_info.chr_start, brca1_info.chr_end, "external_tools.json",write_dir=".cache")
    brca1_clinvar_variants.to_json(".cache/BRCA1_clinvar.json")

def test_queryGnomAD():
    gnomAD_tst_out = queryGnomAD('GRCh38', '17', 43045681, 43124096,"HGNC:1100",'external_tools.json',write_dir=".cache")
    gnomAD_tst_out.to_json(".cache/gnomAD_test.json")

def test_querySpliceAI():
    spliceAI_tst_out = querySpliceAI('GRCh38', '17', 43045681, 43124096,'external_tools.json',write_dir=".cache")
    spliceAI_tst_out.to_json(".cache/spliceAI_test.json")
