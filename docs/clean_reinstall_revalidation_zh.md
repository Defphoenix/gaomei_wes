# WES 全新安装与真实配对样本复测手册

本文用于从 GitHub 干净版本重新安装环境，并在独立目录复测一组 tumor-normal
样本。建议创建新的版本化环境根目录，不删除正在使用的环境和既有分析结果。

## 1. 获取干净代码

```bash
cd /PUBLIC/gomics/guofenghua/project
git clone git@github.com:Defphoenix/gaomei_wes.git gaomei_wes_retest
cd gaomei_wes_retest
git log -1 --oneline
bash scripts/run_code_tests.sh
```

已有 Git 仓库时可使用：

```bash
git pull --ff-only origin main
bash scripts/run_code_tests.sh
```

## 2. 创建独立复测环境

不要先删除当前 `/PUBLIC/gomics/guofenghua/envs/wes`。为本次复测创建新 prefix：

```bash
ENV_ROOT=/PUBLIC/gomics/guofenghua/envs/wes_retest_202607

bash scripts/create_conda_envs.sh \
  --env-root "${ENV_ROOT}" \
  --mamba-bin mamba \
  --production
```

安装器将创建核心、VEP、SnpEff、MHCflurry、HLA*LA 和 CNVkit 环境，保存
manifest，并自动执行轻量回归。Manta 不属于默认生产安装。

检查关键工具：

```bash
source "${ENV_ROOT}/env.sh"
gatk --version
samtools --version
mosdepth --version
msisensor-pro 2>&1 | head
bash scripts/run_vep_env.sh --help >/dev/null
bash scripts/run_snpeff_env.sh -version
bash scripts/run_cnvkit_env.sh version
```

## 3. 验证独立资源

软件环境不包含参考数据库。至少确认：

```bash
REFERENCE_DIR=/path/to/reference_data
REFERENCE_FASTA=${REFERENCE_DIR}/hg38/Homo_sapiens_assembly38.fasta
CAPTURE_BED=/path/to/capture_targets.bed

test -s "${REFERENCE_FASTA}"
test -s "${REFERENCE_FASTA}.fai"
test -s "${REFERENCE_FASTA%.fasta}.dict"
test -s "${REFERENCE_DIR}/dbsnp_146.hg38.vcf.gz"
test -s "${REFERENCE_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
test -d "${REFERENCE_DIR}/vep_cache/homo_sapiens/115_GRCh38"
test -s "${REFERENCE_DIR}/protein/protein.fa"
test -s "${CAPTURE_BED}"
```

HLA*LA graph、MHCflurry models、MSIsensor list、Mutect2 common SNP、gnomAD
AF-only 和 PoN 需要分别检查。缺少后面三项时流程可以研发运行，但会在配置检查中
明确警告，不能标记为完整生产级体细胞过滤。

## 4. 创建独立复测项目

使用软链接，不复制 FASTQ：

```bash
OUT_DIR=/PUBLIC/gomics/guofenghua/project/wes_revalidation/pair01_retest

bash scripts/create_wes_project.sh \
  --mode tumor-normal \
  --tumor-fastq-source /path/to/tumor_fastq_dir \
  --normal-fastq-source /path/to/normal_fastq_dir \
  --out-dir "${OUT_DIR}" \
  --tumor-id tumor01 \
  --normal-id normal01 \
  --project-name pair01_retest \
  --copy-mode link \
  --reference-dir "${REFERENCE_DIR}" \
  --reference-genome "${REFERENCE_FASTA}" \
  --interval-bed "${CAPTURE_BED}" \
  --env-root "${ENV_ROOT}" \
  --no-testdata
```

运行前检查实际生效配置：

```bash
cd "${OUT_DIR}"
bash run_pipeline.sh check
grep -E 'TUMOR_SAMPLE_ID=|NORMAL_SAMPLE_ID=|TUMOR_BAM=|NORMAL_BAM=|CNV_METHOD=' \
  configs/*.somatic.config.sh
```

## 5. 运行与监控

前台运行：

```bash
bash run_pipeline.sh
```

稳定 SSH 会话或调度系统中运行更合适。实时监控：

```bash
tail -f logs/pipeline.log
```

单步恢复仍然保留：

```bash
bash run_pipeline.sh status
bash run_pipeline.sh from somatic 7c
bash run_pipeline.sh step somatic 8
```

## 6. 验收门槛

复测完成后至少检查：

1. Tumor/normal BQSR BAM 均通过 `samtools quickcheck`。
2. Mutect2 VCF header 的 `tumor_sample`、`normal_sample` 与配置一致。
3. `mutect2.filtered.vcf.gz` 与 `mutect2.pass.vcf.gz` 可被 bcftools 读取。
4. F1R2 orientation model 成功；有 common SNP 时 contamination table 成功。
5. VEP VCF、neoantigen FASTA/TSV、HLA binding 表均存在且非空。
6. MSI score 必须位于 0-100%，并与 unstable/total 位点比例一致。
7. CNVkit 只有在 matched/pooled normal 基线存在时才生成 CNR/CNS；否则只能有
   `depth_qc.tsv`，不得报告为正式 CNV。
8. TMB 报告同时保留 accepted/rejected 表和分母来源。
9. 最终报告不得把研发阈值结果描述为临床结论。

完成单对样本复测后，再进入多 normal 的 PoN/CNV reference 构建，以及不同混合
比例、不同深度和不同 VAF 的 LOD 实验。
