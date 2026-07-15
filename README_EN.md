# WES Tumor-Normal and Single-Sample Pipeline

[中文](README.md) | **English**

This repository provides a server-deployable Shell/Python WES workflow with
single-sample germline and matched tumor-normal somatic modes. It is suitable
for development, benchmarking, and workflow validation. Clinical or industrial
deployment still requires cohort validation, calibrated reference resources,
and a governed interpretation layer.

## Workflow

```text
FASTQ -> FastQC -> fastp -> BWA-MEM -> sort/index -> duplicate marking
      -> post-alignment QC -> optional BQSR
      -> HaplotypeCaller or matched Mutect2 -> filtering
      -> optional SnpEff/VEP -> neoantigen peptides -> optional HLA binding
      -> optional CNV/MSI/SV -> coverage/TMB -> summary/MultiQC
```

The project supports:

- Single-sample germline calling with HaplotypeCaller.
- Matched tumor-normal somatic calling with Mutect2.
- One-command project execution plus `step` and `from` debugging.
- VEP 115 offline annotation and 8-15mer neoantigen peptide generation.
- Optional MHCflurry/NetMHCpan binding prediction when HLA alleles are supplied.
- Optional CNVkit, MSIsensor-pro, Manta, coverage, TMB, and summary modules.

It does not yet provide cohort joint genotyping, automatic HLA typing, a full
Mutect2 contamination/orientation-bias workflow, or a validated clinical report.
See [the Chinese audit and roadmap](docs/pipeline_audit_zh.md) for details.

## Install With Mamba

```bash
git clone git@github.com:Defphoenix/gaomei_wes.git
cd gaomei_wes

bash scripts/create_conda_envs.sh \
  --env-root /PUBLIC/gomics/guofenghua/envs/wes \
  --mamba-bin mamba \
  --with-hla \
  --with-cnv
```

Created environments:

| Prefix | Purpose |
|---|---|
| `big_wes_pipeline_env` | Core QC, alignment, GATK, VCF, MSI, and reporting tools |
| `wes_vep_env` | Ensembl VEP 115 |
| `wes_hla_env` | Optional MHCflurry binding prediction |
| `wes_cnv_env` | Optional CNVkit analysis |
| `wes_sv_env` | Optional archived Manta 1.6.0; add `--with-sv` on Linux |

The installer writes the resolved package list and platform-specific explicit
lock files to `ENV_ROOT/manifests/`.

After a future `git pull`, add `--update-existing` to apply changed YML files
to existing prefixes. Without that flag, existing environments are left intact.

Load the environment paths and test the installation:

```bash
source /PUBLIC/gomics/guofenghua/envs/wes/env.sh
gatk --version
vep --help
```

MHCflurry models are downloaded separately:

```bash
mamba run -p /PUBLIC/gomics/guofenghua/envs/wes/wes_hla_env \
  mhcflurry-downloads fetch models_class1_presentation
```

## Reference Layout

```text
reference_data/
  hg38/Homo_sapiens_assembly38.fasta[.fai]
  hg38/Homo_sapiens_assembly38.dict
  dbsnp_146.hg38.vcf.gz[.tbi]
  Mills_and_1000G_gold_standard.indels.hg38.vcf.gz[.tbi]
  1000G_phase1.snps.high_confidence.hg38.vcf.gz[.tbi]
  vep_cache/homo_sapiens/115_GRCh38/
  protein/protein.fa
  msisensor/hg38.list
  capture_targets.bed
```

Use a capture BED matching the actual WES kit. A Mutect2 panel of normals,
population germline resource, CNVkit reference, and MSI baseline are
project/cohort resources rather than universal files.

## Create And Run A Matched Project

```bash
bash scripts/create_wes_project.sh \
  --mode tumor-normal \
  --tumor-fastq-source /data/TUMOR01 \
  --normal-fastq-source /data/NORMAL01 \
  --out-dir /analysis/TUMOR01_vs_NORMAL01 \
  --tumor-id TUMOR01 \
  --normal-id NORMAL01 \
  --copy-mode link \
  --reference-dir /reference_data \
  --reference-genome /reference_data/hg38/Homo_sapiens_assembly38.fasta \
  --interval-bed /reference_data/capture_targets.bed \
  --env-root /PUBLIC/gomics/guofenghua/envs/wes

cd /analysis/TUMOR01_vs_NORMAL01
bash run_pipeline.sh
```

Debug and resume commands:

```bash
bash run_pipeline.sh check
bash run_pipeline.sh status somatic
bash run_pipeline.sh step somatic 7d
bash run_pipeline.sh from somatic 7c
```

For a single germline sample, use `--mode single`, `--fastq-source`, and
`--sample-id`, then run the generated root `run_pipeline.sh` in the same way.

## Important Boundaries

- NetMHCpan, ANNOVAR, and COSMIC require their own registration or license flow
  and are not downloaded automatically.
- `simple` HLA scoring is an explicit smoke-test mode only. `auto` fails when
  neither NetMHCpan nor MHCflurry is available.
- The mosdepth CNV fallback is a development estimate, not a replacement for a
  calibrated CNVkit/FACETS workflow.
- MSI without a compatible site list/baseline is not a formal MSI call.
- MultiQC is a QC aggregation report, not a clinical interpretation report.
