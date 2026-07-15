# WES 流程首次安装与版本更新指南

本文适用于 Linux 分析服务器，推荐使用 `mamba` 和绝对路径 prefix 管理环境。
环境、参考数据和分析项目彼此独立：更新代码不会删除结果，更新环境也不会修改
FASTQ、BAM、VCF 或参考数据库。

## 1. 首次安装

### 1.1 获取代码

```bash
git clone git@github.com:Defphoenix/gaomei_wes.git
cd gaomei_wes
```

没有配置 GitHub SSH 时使用 HTTPS：

```bash
git clone https://github.com/Defphoenix/gaomei_wes.git
cd gaomei_wes
```

如果代码是从本地 Mac 复制到服务器，可使用兼容 macOS 旧版 rsync 的命令：

```bash
rsync -avP \
  --exclude '.git/' \
  --exclude 'results/' \
  --exclude 'logs/' \
  --exclude '.runtime/' \
  /local/path/gaomei_wes/ \
  user@server:/PUBLIC/project/wes_pipeline/
```

不要使用 macOS rsync 2.6.9 不支持的 `--info=progress2`。

### 1.2 创建环境

```bash
cd /PUBLIC/project/wes_pipeline

ENV_ROOT=/PUBLIC/gomics/guofenghua/envs/wes

bash scripts/create_conda_envs.sh \
  --env-root "${ENV_ROOT}" \
  --mamba-bin mamba \
  --with-hla \
  --with-hla-typing \
  --with-cnv
```

推荐配置会创建以下 prefix：

| Prefix | 用途 | 关键运行时 |
|---|---|---|
| `big_wes_pipeline_env` | QC、比对、GATK、VCF、MSI、报告 | Java 17 |
| `wes_vep_env` | Ensembl VEP | Perl/VEP 115 |
| `wes_snpeff_env` | SnpEff、SnpSift | Java 21 |
| `wes_hla_env` | MHCflurry binding prediction | Python |
| `wes_hla_typing_env` | HLA*LA 分型 | Linux |
| `wes_cnv_env` | CNVkit | 独立 Python/R 依赖 |

安装器会生成 `${ENV_ROOT}/env.sh`，并将实际解析的软件版本、build 和下载地址
写入 `${ENV_ROOT}/manifests/`。

### 1.3 加载并验证

```bash
source "${ENV_ROOT}/env.sh"

gatk --version
bwa 2>&1 | head
samtools --version
bcftools --version
picard -h 2>&1 | grep -m1 PicardCommandLine
vep --help >/dev/null && echo VEP_OK
bash scripts/run_snpeff_env.sh -version
bash scripts/run_snpsift_env.sh -h 2>&1 | head
mamba run -p "${ENV_ROOT}/wes_cnv_env" cnvkit.py version
```

这些环境由绝对路径创建，激活时也必须使用完整路径：

```bash
mamba activate /PUBLIC/gomics/guofenghua/envs/wes/big_wes_pipeline_env
```

不能只运行 `mamba activate big_wes_pipeline_env`。

## 2. 更新代码和环境

### 2.1 仅更新流程代码

```bash
cd /PUBLIC/project/gaomei_wes
git pull --ff-only origin main
bash -n run_pipeline.sh config.sh scripts/*.sh
```

如果没有 `.git` 目录，说明当前目录是复制版而不是 Git 仓库。继续使用 rsync 覆盖
代码即可，或者重新 `git clone`；这不影响流程本身运行。

### 2.2 代码新增了一个环境

重新执行首次安装命令。已存在的完整环境会显示 `[SKIP]`，安装器只创建缺失的
prefix，然后重新运行工具验证：

```bash
bash scripts/create_conda_envs.sh \
  --env-root "${ENV_ROOT}" \
  --mamba-bin mamba \
  --with-hla \
  --with-hla-typing \
  --with-cnv
```

例如引入 `wes_snpeff_env` 后，不需要删除或重建 GATK 主环境。

### 2.3 YML 修改了已有环境的软件版本

需要让现有 prefix 与最新版 YML 同步时，增加 `--update-existing`：

```bash
bash scripts/create_conda_envs.sh \
  --env-root "${ENV_ROOT}" \
  --mamba-bin mamba \
  --with-hla \
  --with-hla-typing \
  --with-cnv \
  --update-existing
```

该操作会重新求解依赖并清理 YML 中已移除的软件，可能耗时较长。更新完成后检查
`${ENV_ROOT}/manifests/*.packages.txt`，并重新执行工具验证。

### 2.4 更新后重新生成项目

如果新版本增加了 config 变量、分析步骤或独立环境，推荐用最新版
`scripts/create_wes_project.sh` 重新生成项目。不要直接复用很早版本自动生成的
config；FASTQ 可以继续使用软链接，不必重复复制。

## 3. 安装后仍需准备的数据

Conda 环境只包含软件，不自动包含大型数据库或许可数据。至少检查：

```text
reference_data/
├── hg38/Homo_sapiens_assembly38.fasta
├── dbsnp_146.hg38.vcf.gz
├── Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
├── 1000G_phase1.snps.high_confidence.hg38.vcf.gz
├── vep_cache/homo_sapiens/115_GRCh38/
├── protein/protein.fa
├── mutect2/small_exac_common_3.hg38.vcf.gz
├── hla/PRG_MHC_GRCh38_withIMGT/
└── msi/<capture-kit>.list
```

还需要参考 FASTA 的 BWA、samtools 和 GATK 字典索引，以及与捕获试剂匹配的 BED。
MHCflurry 模型和 HLA*LA graph 需要单独下载；MSI list 应根据参考基因组和捕获
BED 构建或筛选。

首次使用 MHCflurry 时下载模型：

```bash
mamba run -p "${ENV_ROOT}/wes_hla_env" \
  mhcflurry-downloads fetch models_class1_presentation
```

HLA*LA graph 下载完成后准备并校验：

```bash
bash scripts/prepare_hlala_graph.sh \
  --archive /path/to/PRG_MHC_GRCh38_withIMGT.tar.gz \
  --reference-dir /path/to/reference_data \
  --env-prefix "${ENV_ROOT}/wes_hla_typing_env"
```

SnpEff 软件环境不等于 SnpEff 基因组数据库。运行步骤 7b 前，确认 config 中的
`SNPEFF_DB`、`SNPEFF_DATA_DIR` 和 `SNPEFF_CONFIG` 指向服务器上的同一套 GRCh38
数据库；VEP 同样需要与 `VEP_CACHE_VERSION` 对应的离线 cache。

## 4. 常见安装问题

### GATK 提示 class file version 61

GATK 4.6 使用主环境 Java 17。先加载：

```bash
source "${ENV_ROOT}/env.sh"
java -version
gatk --version
```

### SnpEff 提示 class file version 65

新版 SnpEff 需要 Java 21。不要升级 GATK 主环境的 Java；使用独立启动器：

```bash
bash scripts/run_snpeff_env.sh -version
```

### BWA 已存在但旧安装器显示 FAILED

旧验证器使用登录 shell，可能重置 `PATH`。最新版已经修复。手动确认：

```bash
"${ENV_ROOT}/big_wes_pipeline_env/bin/bwa" 2>&1 | head
```

### Picard 显示帮助后被判定 FAILED

`picard -h` 能正常列出程序但可能返回非零状态。最新版验证器检查帮助内容，
不再只判断退出码：

```bash
picard -h 2>&1 | grep -m1 PicardCommandLine
```

### 环境创建中断

先用 `mamba env list` 判断 prefix 是否完整。安装器会识别包含
`conda-meta/history` 的完整环境；不完整目录应在确认没有需要保留的数据后删除，
再重新执行安装命令。

## 5. 推荐更新顺序

```text
备份项目 config
→ git pull 或 rsync 同步代码
→ bash -n 静态检查
→ 运行 create_conda_envs.sh
→ 检查 manifests 和工具版本
→ 检查参考数据路径
→ 重新生成测试项目
→ 小样本端到端测试
→ 正式样本运行
```

更完整的项目生成、混样和服务器运行示例见
[`server_deployment_and_mix_test.md`](server_deployment_and_mix_test.md)。
