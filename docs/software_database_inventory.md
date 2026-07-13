# WES Pipeline Software and Database Inventory

Generated for local workspace: `/Users/mac/Documents/wes/test_run_workspace`

## Local Status Summary

| Category | Item | Current local status | Suggested action |
|---|---|---:|---|
| Core reference | GRCh38 FASTA + `.fai` + BWA index + `.dict` | Installed at `/Users/mac/Documents/wes/reference_data/hg38` | Keep |
| Germline known sites | dbSNP 146 hg38 + index | Installed | Keep |
| Germline known indels | Mills and 1000G gold standard indels + index | Installed | Keep |
| Germline known SNPs | 1000G phase1 high-confidence SNPs + index | Installed | Keep |
| VEP cache | homo_sapiens VEP 115 GRCh38 cache | Extracted at `/Users/mac/Documents/wes/reference_data/vep_cache/homo_sapiens/115_GRCh38` | Keep |
| VEP executable | `vep` | Installed at `/Users/mac/Documents/wes/.conda_envs/wes_vep_env/bin/vep` via wrapper `scripts/run_vep_env.sh` | Keep; tested with step 7c |
| ANNOVAR software | `table_annovar.pl`, `annotate_variation.pl` | Missing | Manual licensed download |
| ANNOVAR databases | `humandb/hg38_*` | Missing | Download after ANNOVAR installed |
| SnpEff executable | `snpEff` / `SnpSift` | Installed in `/Users/mac/anaconda3/envs/big_wes_pipeline_env` | Keep |
| SnpEff database | `GRCh38.105` database | Not found under `/Users/mac/Documents/wes/reference_data` | Download only if SnpEff path is used |
| CNVkit | `cnvkit.py` | Missing | Optional; pipeline has mosdepth fallback |
| MSI | `msisensor-pro` | Installed | Provide list/baseline for formal MSI |
| HLA binding | `netMHCpan`, `mhcflurry-predict` | Missing | netMHCpan manual license or install mhcflurry |
| Neoantigen protein FASTA | Ensembl/GENCODE protein FASTA | Missing | Download if neoantigen module is used |

## Software Table

| Priority | Tool | Purpose | Current status | Install route | Notes |
|---:|---|---|---:|---|---|
| P0 | FastQC | raw read QC | Installed | conda/bioconda | Core QC |
| P0 | fastp | trimming | Installed | conda/bioconda | Used by step 2 |
| P0 | BWA | alignment | Installed | conda/bioconda | Used by step 3 |
| P0 | Samtools | BAM/FASTA operations | Installed | conda/bioconda | Core |
| P0 | GATK4 | variant calling/BQSR/filtering | Installed | conda/bioconda | Core |
| P0 | BCFtools/HTSlib | VCF operations/indexing | Installed | conda/bioconda | Core |
| P0 | BEDtools | BED interval operations | Installed | conda/bioconda | Core |
| P0 | Picard | duplicate marking/QC | Installed | conda/bioconda | Core |
| P1 | VEP 115 | functional annotation | Installed and tested | `wes_vep_env.yml` | Cache already present |
| P1 | ANNOVAR | functional/population/clinical annotation | Missing | Manual official download | Requires registration/license |
| P1 | SnpEff/SnpSift | functional annotation | Installed executable | conda/bioconda | DB still needed |
| P1 | MultiQC | QC aggregation | Installed | pip/conda | Already in env |
| P1 | mosdepth | coverage and CNV fallback | Installed | conda/bioconda | Used by CNV fallback |
| P1 | MSIsensor-pro | MSI | Installed | conda/bioconda | Formal run needs microsatellite list/baseline |
| P2 | CNVkit | CNV | Missing | conda/bioconda | Optional; preferred for formal CNV |
| P2 | Manta | SV | Missing | conda/bioconda or separate env | Optional SV module |
| P2 | netMHCpan | HLA binding | Missing | Manual DTU licensed download | Best known option, license required |
| P2 | mhcflurry | HLA binding | Missing | `wes_hla_env.yml` | Open install route, models may require extra downloads |
| P2 | OptiType / arcasHLA | HLA typing | Missing | separate env | Needed if HLA alleles are not provided |

## Database Table

| Priority | Database/resource | Module | Current status | Approx size | Download / preparation route | Notes |
|---:|---|---|---:|---:|---|---|
| P0 | GRCh38 reference FASTA | alignment/calling | Installed | ~3-4 GB plus indexes | Already local | Keep under `reference_data/hg38` |
| P0 | dbSNP hg38 VCF | BQSR/annotation | Installed | several GB | Already local | `dbsnp_146.hg38.vcf.gz` |
| P0 | Mills indels hg38 VCF | BQSR | Installed | small-medium | Already local | Used by BQSR |
| P0 | 1000G high-confidence SNPs | BQSR | Installed | medium | Already local | Optional known site |
| P1 | VEP cache 115 GRCh38 | VEP | Installed/extracted | 24 GB extracted | Already local | Cache metadata includes COSMIC 101, dbSNP 156, ClinVar 202502, 1000G phase3, gnomAD v4.1 |
| P1 | VEP plugins data: dbNSFP | VEP | Missing | large, often 30-100+ GB | manual/plugin-specific | Supports SIFT/PolyPhen-like aggregated scores, REVEL etc. |
| P1 | VEP plugins data: CADD | VEP | Missing | large | CADD release files | Optional, large |
| P1 | VEP plugins data: REVEL | VEP | Missing | medium-large | REVEL release | Optional |
| P1 | ANNOVAR humandb `refGene` | ANNOVAR | Missing | small | `annotate_variation.pl -buildver hg38 -downdb refGene humandb/` | Basic gene annotation |
| P1 | ANNOVAR humandb `cytoBand` | ANNOVAR | Missing | small | ANNOVAR download | Cytoband annotation |
| P1 | ANNOVAR humandb `avsnp150/151` | ANNOVAR | Missing | medium | ANNOVAR download | dbSNP ID annotation |
| P1 | ANNOVAR humandb `clinvar_*` | ANNOVAR | Missing | medium | ANNOVAR download | Clinical annotation |
| P1 | ANNOVAR humandb `gnomad*` | ANNOVAR | Missing | large | ANNOVAR download | Population AF filtering |
| P1 | ANNOVAR humandb `dbnsfp*` | ANNOVAR | Missing | very large | ANNOVAR download | Functional prediction scores |
| P1 | SnpEff GRCh38 database | SnpEff | Missing | several GB | `snpEff download GRCh38.105` | Matches current `SNPEFF_DB=GRCh38.105`; alternatives include `GRCh38.p14`, MANE, and UCSC `hg38` |
| P1 | ClinVar VCF hg38 | SnpSift/bcftools annotate | Missing standalone | medium | NCBI ClinVar VCF | Useful even without ANNOVAR |
| P1 | COSMIC coding/noncoding | cancer annotation | Missing | licensed | COSMIC account/license | Cannot auto-download without credentials |
| P1 | MSIsensor-pro microsatellite list | MSI | Missing | small-medium | `msisensor-pro scan` or prebuilt list | Formal MSI needs this |
| P1 | MSIsensor-pro baseline | MSI tumor-only | Missing | medium | build/download baseline | Tumor-only formal calling improves with baseline |
| P1 | Ensembl/GENCODE protein FASTA | neoantigen | Missing | small-medium | Ensembl/GENCODE | Needed to reconstruct peptide context |
| P2 | Panel of Normals | Mutect2/CNV/MSI | Missing | project-specific | build from normals | Strongly recommended for somatic WES |
| P2 | CNVkit pooled reference `.cnn` | CNVkit | Missing | project-specific | build from normals | Formal CNV improves with matched normals |
| P2 | target capture BED | WES coverage/CNV | Using test BED only | small | vendor kit BED | Required for real WES |

## Recommended Local Directory Layout

```text
/Users/mac/Documents/wes/reference_data/
  hg38/
    Homo_sapiens_assembly38.fasta
    Homo_sapiens_assembly38.fasta.fai
    Homo_sapiens_assembly38.dict
  known_sites/
    dbsnp_146.hg38.vcf.gz
    Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    1000G_phase1.snps.high_confidence.hg38.vcf.gz
  vep_cache/
    homo_sapiens/115_GRCh38/
  annovar/
    annovar/
    humandb/
  snpeff/
    data/
  clinvar/
  cosmic/
  msisensor/
  protein/
  hla/
```

## Suggested Download Order

1. Download SnpEff GRCh38 database if we keep SnpEff as a parallel annotation path.
2. Download/install ANNOVAR manually, then populate `annovar/humandb`.
3. Add ClinVar/gnomAD/dbNSFP either through ANNOVAR humandb or VEP plugins; avoid duplicating huge databases unless needed.
4. Prepare MSIsensor-pro list/baseline.
5. Add HLA binding tool: netMHCpan if licensed, otherwise mhcflurry.
6. Optional: CNVkit and Manta in separate optional environments.

## Licensing / Manual-download Items

| Item | Reason |
|---|---|
| ANNOVAR | Official download requires registration and license acceptance. |
| COSMIC | Requires COSMIC account/license. |
| netMHCpan | Requires DTU academic/commercial license flow. |
| Some CADD/dbNSFP/REVEL resources | Large files and may require accepting source-specific terms. |
