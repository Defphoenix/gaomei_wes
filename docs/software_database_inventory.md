# WES 软件与数据库清单

本清单描述 GitHub 版本的依赖边界，不记录某一台电脑的绝对路径。安装后真实版本
保存在 `ENV_ROOT/manifests/*.packages.txt`，可用于服务器复现和审计。

## 自动安装的软件

| 环境 | 软件 | 用途 |
|---|---|---|
| 核心 | Python 3.10、OpenJDK 17 | 脚本与 GATK 运行时 |
| 核心 | FastQC、fastp、Seqtk、Trimmomatic | FASTQ QC、修剪和抽样 |
| 核心 | BWA、Samtools、BCFtools、HTSlib、BEDtools | 比对及 BAM/VCF/BED 操作 |
| 核心 | GATK4、Picard | 胚系/体细胞检测、BQSR、重复标记 |
| 核心 | SnpEff/SnpSift | 可选功能注释；数据库另行准备 |
| 核心 | mosdepth、MSIsensor-pro、MultiQC | 深度、MSI 与 QC 汇总 |
| VEP | Ensembl VEP 115.2 | 离线转录本和蛋白后果注释 |
| HLA，可选 | MHCflurry | HLA-I binding；模型另行下载 |
| HLA typing，可选、Linux | HLA*LA 1.0.4 | WES BAM G-group分型；graph另行准备 |
| CNV，可选 | CNVkit 0.9.12 | 正式 CNV 主工具 |
| SV，可选、Linux | Manta 1.6.0 | 兼容 SV 模块；上游项目已归档 |

推荐安装：

```bash
bash scripts/create_conda_envs.sh \
  --env-root /PUBLIC/gomics/guofenghua/envs/wes \
  --mamba-bin mamba \
  --with-hla \
  --with-hla-typing \
  --with-cnv
```

Linux 服务器确实需要 Manta 时增加 `--with-sv`。首次安装 MHCflurry 模型可增加
`--fetch-hla-models`，但该步骤会联网并下载额外数据。

## 不能自动下载的软件或数据

| 项目 | 原因 | 流程中的位置 |
|---|---|---|
| NetMHCpan | DTU 许可和申请流程 | 可替代 MHCflurry 做 binding |
| ANNOVAR 与 humandb | 注册/许可要求 | 可选并行注释；新抗原表可合并其结果 |
| COSMIC | 账户和许可要求 | 肿瘤知识注释 |
| 部分 CADD/dbNSFP/REVEL 数据 | 许可、版本和体积差异 | VEP 插件 |

VEP cache 本身不等于自动拥有所有插件数据；`--everything` 也不会凭空提供
CADD/REVEL。需要这些评分时，应明确下载对应版本并配置插件参数。

## 参考数据

| 优先级 | 数据 | 用途 | 说明 |
|---:|---|---|---|
| P0 | GRCh38 FASTA、FAI、DICT、BWA index | 全流程 | contig 命名必须与 BED/VCF 一致 |
| P0 | 捕获试剂盒 BED | WES 调用、覆盖度、CNV、TMB | 使用厂商对应版本，不使用通用测试 BED |
| P0 | VEP 115 GRCh38 cache | VEP | 目录应为 `homo_sapiens/115_GRCh38` |
| P0 | Ensembl 115 `pep.all` FASTA | 新抗原 | 蛋白 ID 需匹配 VEP ENSP/Feature |
| P1 | dbSNP、Mills、1000G known-sites | BQSR | 全部使用同一 GRCh38 build |
| P1 | gnomAD AF-only resource | Mutect2 | 用于germline prior |
| P1 | common biallelic SNP VCF | Mutect2 | `GetPileupSummaries`/污染估计，需索引 |
| P1 | Panel of Normals | Mutect2 | 同平台、同试剂盒 normal 构建 |
| P1 | CNVkit normal reference | CNV | 同批次 normal 建立并版本化 |
| P1 | MSIsensor site list/baseline | MSI | 应按参考版本、BED 和平台校准 |
| P2 | ClinVar VCF | 临床注释 | 可通过 VEP custom 或其他注释器接入 |
| P2 | SnpEff GRCh38 database | SnpEff | 只有启用 SnpEff 时需要 |
| P1 | HLA*LA PRG graph | HLA分型 | 约2.3GB压缩包，prepare时约需40GB内存 |
| P1 | effective coding BED | TMB分母 | capture BED与CDS求交并合并，按panel版本化 |

## 推荐目录

```text
reference_data/
  hg38/
    Homo_sapiens_assembly38.fasta
    Homo_sapiens_assembly38.fasta.fai
    Homo_sapiens_assembly38.dict
  known_sites/
  mutect2/
    gnomad.af-only.vcf.gz
    small_exac_common_3.hg38.vcf.gz
    panel_of_normals.vcf.gz
  hla/PRG_MHC_GRCh38_withIMGT/
  tmb/effective_coding_regions.bed
  vep_cache/homo_sapiens/115_GRCh38/
  protein/protein.fa
  cnvkit/reference.cnn
  msisensor/hg38.capture.list
  capture_targets.bed
```

## 重要说明

- MHCflurry/NetMHCpan 做的是 binding prediction，不等于 HLA typing。
- mosdepth 深度比例是 CNV 连通性 fallback，不应替代校准后的 CNVkit 结果。
- 缺少 MSI 位点表或基线时只能产生测试状态，不能作为正式 MSI 结论。
- 数据库版本、下载日期、校验值和许可信息都应进入项目 provenance。
