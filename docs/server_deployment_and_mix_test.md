# WES Pipeline Server Deployment and Mixed-sample Test

This document describes how to deploy the GitHub version on a server, create
mamba/conda environments, generate mixed FASTQ benchmark data, and create
runnable sample configs.

## 1. Get Or Update The Code

Recommended server workflow:

```bash
mkdir -p /PUBLIC/gomics/guofenghua/project
cd /PUBLIC/gomics/guofenghua/project

git clone git@github.com:Defphoenix/gaomei_wes.git wes_pipeline
cd wes_pipeline
```

If the server does not have GitHub SSH configured, use HTTPS instead:

```bash
git clone https://github.com/Defphoenix/gaomei_wes.git wes_pipeline
cd wes_pipeline
```

After local changes are pushed to GitHub, update the server copy with:

```bash
cd /PUBLIC/gomics/guofenghua/project/wes_pipeline
git pull
bash scripts/create_conda_envs.sh \
  --env-root /PUBLIC/gomics/guofenghua/envs/wes \
  --mamba-bin mamba \
  --with-hla \
  --with-cnv \
  --update-existing
```

If you cannot use GitHub on the server, rsync is still possible, but exclude
runtime outputs:

```bash
rsync -a \
  --exclude results/ \
  --exclude logs/ \
  --exclude cache/ \
  --exclude .runtime/ \
  /Users/mac/Documents/wes/test_run_workspace/ \
  user@server:/PUBLIC/gomics/guofenghua/project/wes_pipeline/
```

The repository contains one tiny paired FASTQ for display/smoke tests:

```text
testdata/variant_sim/sim100_R1.fastq.gz
testdata/variant_sim/sim100_R2.fastq.gz
```

## 2. Create Environments

On the server:

```bash
cd /PUBLIC/gomics/guofenghua/project/wes_pipeline

# Prefer mamba. Use --mamba-bin micromamba if that is what the server provides.
bash scripts/create_conda_envs.sh \
  --env-root /PUBLIC/gomics/guofenghua/envs/wes \
  --mamba-bin mamba \
  --with-hla \
  --with-cnv
```

If your server has global channels such as `r` in `.condarc`, keep the default
`--override-channels --channel-priority flexible` behavior. It prevents unrelated
global channels from entering the solve.

Picard is installed in the core environment. CNVkit is isolated in
`wes_cnv_env` because its R/Python dependencies are much larger. Manta is also
isolated and is only created with `--with-sv`; use it on Linux only. The
installer writes exact resolved versions under `ENV_ROOT/manifests/`.

If a previous failed environment directory exists, either remove it manually
after confirming the path, or let the installer clean only incomplete prefix
directories:

```bash
bash scripts/create_conda_envs.sh \
  --env-root /PUBLIC/gomics/guofenghua/envs/wes \
  --mamba-bin mamba \
  --with-hla \
  --clean-incomplete
```

Then use these paths in generated configs:

```bash
MAIN_ENV_PREFIX=/PUBLIC/gomics/guofenghua/envs/wes/big_wes_pipeline_env
VEP_ENV_PREFIX=/PUBLIC/gomics/guofenghua/envs/wes/wes_vep_env
HLA_ENV_PREFIX=/PUBLIC/gomics/guofenghua/envs/wes/wes_hla_env
CNV_ENV_PREFIX=/PUBLIC/gomics/guofenghua/envs/wes/wes_cnv_env
SV_ENV_PREFIX=/PUBLIC/gomics/guofenghua/envs/wes/wes_sv_env
```

The optional HLA environment installs MHCflurry through conda/bioconda. It does
not install NetMHCpan. NetMHCpan is a DTU standalone tool and should be manually
downloaded/installed only after accepting the appropriate license. In configs,
`HLA_BINDING_TOOL=auto` tries NetMHCpan first and then MHCflurry. If neither is
available it stops; the built-in `simple` mode must be selected explicitly and
is only a smoke test.

MHCflurry separates the software package and prediction models. After creating
the HLA environment, fetch the class-I presentation models once:

```bash
mamba run -p /PUBLIC/gomics/guofenghua/envs/wes/wes_hla_env \
  mhcflurry-downloads fetch models_class1_presentation
```

These environments are prefix-based, so activate by full path, not by short
name:

```bash
mamba activate /PUBLIC/gomics/guofenghua/envs/wes/big_wes_pipeline_env
```

For command-line testing, the easiest method is to source the generated helper:

```bash
source /PUBLIC/gomics/guofenghua/envs/wes/env.sh
which java
java -version
gatk --version
```

GATK 4.6 requires Java 17. If you test GATK manually from a shell, put the main
environment first in PATH and set JAVA_HOME:

```bash
export MAIN_ENV_PREFIX=/PUBLIC/gomics/guofenghua/envs/wes/big_wes_pipeline_env
export JAVA_HOME=${MAIN_ENV_PREFIX}
export PATH=${MAIN_ENV_PREFIX}/bin:$PATH
java -version
gatk --version
```

## 3. Generate Mixed FASTQ Benchmark Data

This replaces the earlier hard-coded HG002/HG003/HG004 script. It keeps the same
logic but makes sample directories, ratios, pairs, seed, and output root explicit.

```bash
cd /PUBLIC/gomics/guofenghua/project/wes_pipeline

export PATH=/PUBLIC/gomics/guofenghua/envs/wes/big_wes_pipeline_env/bin:$PATH

bash scripts/make_mix_fastq.sh \
  --base-dir /PUBLIC/gomics/guofenghua/project/project_01_wespipline/hg002_benmark \
  --out-base /PUBLIC/gomics/guofenghua/project/project_01_wespipline/hg002_benmark/moni \
  --sample HG002:HG002_Sample_2A1 \
  --sample HG003:HG003_Sample_3A1 \
  --sample HG004:HG004_Sample_4A1 \
  --pairs HG002:HG003,HG002:HG004 \
  --ratios 5,10,15,20,30 \
  --seed 100
```

For a small first test, generate only a fraction of the full mixed data. The
example below writes to `moni_small` and uses `--total-fraction 0.01`, which is
usually closer to a few hundred MB per mixed FASTQ set than the full benchmark.
Use `0.1` for one tenth of the full mix after this small test is validated.

```bash
cd /PUBLIC/gomics/guofenghua/project/wes_pipeline
source /PUBLIC/gomics/guofenghua/envs/wes/env.sh

bash scripts/make_mix_fastq.sh \
  --base-dir /PUBLIC/gomics/guofenghua/project/project_01_wespipline/hg002_benmark \
  --out-base /PUBLIC/gomics/guofenghua/project/project_01_wespipline/hg002_benmark/moni_small \
  --sample HG002:HG002_Sample_2A1 \
  --sample HG003:HG003_Sample_3A1 \
  --sample HG004:HG004_Sample_4A1 \
  --pairs HG002:HG003,HG002:HG004 \
  --ratios 5,10,15,20,30 \
  --total-fraction 0.01 \
  --seed 100

du -sh /PUBLIC/gomics/guofenghua/project/project_01_wespipline/hg002_benmark/moni_small/*
cat /PUBLIC/gomics/guofenghua/project/project_01_wespipline/hg002_benmark/moni_small/mixed_fastq_manifest.tsv
```

Outputs:

```text
moni/mix_HG002_HG003_5pct/mix_R1.fastq.gz
moni/mix_HG002_HG003_5pct/mix_R2.fastq.gz
moni/mixed_fastq_manifest.tsv
```

## 4. Generate A Runnable Single-sample Project

Example using one mixed sample:

```bash
cd /PUBLIC/gomics/guofenghua/project/wes_pipeline

bash scripts/create_wes_project.sh \
  --mode single \
  --fastq-source /PUBLIC/gomics/guofenghua/project/project_01_wespipline/hg002_benmark/moni/mix_HG002_HG003_5pct \
  --out-dir /PUBLIC/gomics/guofenghua/project/wes_runs/mix_HG002_HG003_5pct \
  --sample-id mix_HG002_HG003_5pct \
  --copy-mode link \
  --reference-dir /PUBLIC/gomics/guofenghua/reference_data \
  --reference-genome /PUBLIC/gomics/guofenghua/reference_data/hg38/Homo_sapiens_assembly38.fasta \
  --interval-bed /PUBLIC/gomics/guofenghua/reference_data/capture_targets.bed \
  --env-root /PUBLIC/gomics/guofenghua/envs/wes
```

The generated project contains:

```text
mix_HG002_HG003_5pct/
  pipeline/                 copied pipeline code
  configs/sample.config.sh  generated config
  data/sample/              copied or linked FASTQ
  results/sample/           output root
  run_sample.sh             one-command runner
```

Run it through the stable project-root entry point:

```bash
cd /PUBLIC/gomics/guofenghua/project/wes_runs/mix_HG002_HG003_5pct
bash run_pipeline.sh
```

Single-step debug:

```bash
cd /PUBLIC/gomics/guofenghua/project/wes_runs/mix_HG002_HG003_5pct/pipeline
bash run_pipeline.sh --config ../configs/mix_HG002_HG003_5pct.config.sh check
bash run_pipeline.sh --config ../configs/mix_HG002_HG003_5pct.config.sh step 7c
bash run_pipeline.sh --config ../configs/mix_HG002_HG003_5pct.config.sh status
```

## 4a. Generate A Tumor-normal Paired Project

For standard tumor WES, use tumor-normal mode. This aligns tumor and normal
separately, then runs Mutect2 with the matched normal BAM and filters somatic
variants with `FilterMutectCalls`.

```bash
cd /PUBLIC/gomics/guofenghua/project/wes_pipeline

bash scripts/create_wes_project.sh \
  --mode tumor-normal \
  --tumor-fastq-source /PUBLIC/gomics/guofenghua/project/sample_fastq/TUMOR01 \
  --normal-fastq-source /PUBLIC/gomics/guofenghua/project/sample_fastq/NORMAL01 \
  --out-dir /PUBLIC/gomics/guofenghua/project/wes_runs/TUMOR01_vs_NORMAL01 \
  --tumor-id TUMOR01 \
  --normal-id NORMAL01 \
  --copy-mode link \
  --reference-dir /PUBLIC/gomics/guofenghua/reference_data \
  --reference-genome /PUBLIC/gomics/guofenghua/reference_data/hg38/Homo_sapiens_assembly38.fasta \
  --interval-bed /PUBLIC/gomics/guofenghua/reference_data/capture_targets.bed \
  --env-root /PUBLIC/gomics/guofenghua/envs/wes
```

Run the paired project through the stable project-root entry point:

```bash
cd /PUBLIC/gomics/guofenghua/project/wes_runs/TUMOR01_vs_NORMAL01
bash run_pipeline.sh
```

Important outputs:

```text
results/TUMOR01_vs_NORMAL01/variants/TUMOR01_vs_NORMAL01.mutect2.raw.vcf.gz
results/TUMOR01_vs_NORMAL01/variants/TUMOR01_vs_NORMAL01.mutect2.filtered.vcf.gz
results/TUMOR01_vs_NORMAL01/variants/TUMOR01_vs_NORMAL01.mutect2.pass.vcf.gz
results/TUMOR01_vs_NORMAL01/annotation/TUMOR01_vs_NORMAL01.vep.vcf.gz
```

## 4b. Generate Runnable Projects For All Mixed Samples

After `make_mix_fastq.sh` finishes, use its manifest to generate one runnable
project per mixed sample:

```bash
cd /PUBLIC/gomics/guofenghua/project/wes_pipeline

bash scripts/create_projects_from_manifest.sh \
  --manifest /PUBLIC/gomics/guofenghua/project/project_01_wespipline/hg002_benmark/moni/mixed_fastq_manifest.tsv \
  --out-base /PUBLIC/gomics/guofenghua/project/wes_runs/mix_benchmark \
  --copy-mode link \
  --reference-dir /PUBLIC/gomics/guofenghua/reference_data \
  --reference-genome /PUBLIC/gomics/guofenghua/reference_data/hg38/Homo_sapiens_assembly38.fasta \
  --interval-bed /PUBLIC/gomics/guofenghua/reference_data/capture_targets.bed \
  --env-root /PUBLIC/gomics/guofenghua/envs/wes
```

Then run any generated project:

```bash
bash /PUBLIC/gomics/guofenghua/project/wes_runs/mix_benchmark/mix_HG002_HG003_5pct/run_mix_HG002_HG003_5pct.sh
```

Mixed HG002/HG003/HG004 projects are benchmark/simulation projects. They are
single-sample analyses by default and are useful for testing allele fraction or
purity-like behavior. They are not the same as clinical tumor-normal paired WES.

## 5. Smoke Test With Bundled Demo FASTQ

This checks that the project generator works before using large server FASTQ:

```bash
cd /PUBLIC/gomics/guofenghua/project/wes_pipeline

bash scripts/create_wes_project.sh \
  --fastq-source testdata/variant_sim \
  --out-dir /PUBLIC/gomics/guofenghua/project/wes_runs/demo_sim100 \
  --sample-id sim100 \
  --copy-mode copy \
  --reference-dir /PUBLIC/gomics/guofenghua/reference_data \
  --reference-genome /PUBLIC/gomics/guofenghua/reference_data/hg38/Homo_sapiens_assembly38.fasta \
  --interval-bed /PUBLIC/gomics/guofenghua/reference_data/demo_targets.bed \
  --env-root /PUBLIC/gomics/guofenghua/envs/wes
```

Then:

```bash
bash /PUBLIC/gomics/guofenghua/project/wes_runs/demo_sim100/run_sim100.sh
```

## 6. Notes

- `--copy-mode link` is recommended for large benchmark FASTQ to avoid duplicating
  hundreds of GB.
- `--copy-mode copy` is safer for small demo data and fully portable test cases.
- ANNOVAR, COSMIC, and netMHCpan still require manual license/download steps.
- VEP can run from `wes_vep_env`; its cache path must point to the server-side
  `reference_data/vep_cache`.
