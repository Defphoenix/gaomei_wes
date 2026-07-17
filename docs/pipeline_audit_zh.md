# WES 流程审计与更新路线

更新日期：2026-07-16

## 结论

当前仓库已经形成一条可安装、可创建项目、可单步调试的 WES 研发流程。
从 FASTQ 到 BAM、HaplotypeCaller/Mutect2、VEP 和新抗原候选肽的主链可用于
benchmark 与演示。它还不是可以不经验证直接用于临床交付的“工业级流程”。

这次审计修复了以下发布阻断问题：

1. 新生成项目的一键入口改为执行配置启用的完整流程，不再停在 `7c`。
2. `check` 由完整流程自动执行，同时保留 `step` 和 `from` 调试方式。
3. fastp 和 BWA 增加已有结果检查，减少失败重启时的重复计算。
4. VEP 删除不存在的 `--clinvar` 参数；ClinVar 信息来自 cache 的
   `--check_existing` 或显式 `--custom` 文件。
5. HLA binding 的 `auto` 模式找不到正式预测器时跳过；`true` 模式严格失败，不输出模拟分数。
6. 核心、VEP、HLA binding、HLA typing、CNVkit、Manta 环境分离，并保存精确版本manifest。
7. Mutect2污染/方向偏倚模型、HLA*LA分型和严格VEP TMB已接入主流程。
8. 正式项目默认启用BQSR，somatic配置优先使用BQSR BAM并回退到dedup BAM。
9. `check` 新增配置一致性和关键资源检查；项目生成配置显式记录参考根目录、known-sites、Mutect2资源和MSI list。

## 当前模式边界

| 模式 | 当前实现 | 状态 |
|---|---|---|
| 单样本胚系 | FASTQ 到 HaplotypeCaller 与硬过滤 | 可测试；缺少 gVCF 队列联合分型 |
| 单肿瘤样本 | 可通过 Mutect2 tumor-only 配置实现 | 可测试；强依赖 PoN 和 germline resource |
| 肿瘤-正常配对 | 两个样本独立预处理后进入配对 Mutect2 | 当前推荐入口 |
| 多个独立样本 | manifest 批量创建多个项目 | 已支持任务生成，不等于队列联合分析 |
| 胚系多样本联合分型 | GenomicsDBImport/GenotypeGVCFs | 尚未实现 |
| 同一患者多肿瘤/纵向样本 | 多肿瘤联合体细胞模型 | 尚未实现 |

## 模块成熟度

| 模块 | 当前能力 | 生产前必须补充 |
|---|---|---|
| FASTQ QC/修剪 | FastQC + fastp | 建立批次 QC 阈值和失败规则 |
| 比对/去重 | BWA-MEM + Picard，Samtools fallback | 建议升级评估 BWA-MEM2；保留 read group/LIMS 信息 |
| BQSR | GATK BQSR 脚本已存在 | 正式项目默认启用并校验 known-sites |
| 胚系 SNV/InDel | HaplotypeCaller + 硬过滤 | 增加 gVCF、联合分型、VQSR/机器学习过滤策略 |
| 体细胞 SNV/InDel | 配对Mutect2 + 污染估计 + 方向偏倚过滤 | 补充同平台PoN、资源严格preflight、散射并行和真值验证 |
| VEP | VEP 115 离线 cache | 固定 MANE/canonical 选转录本规则；按需接 ClinVar/dbNSFP |
| ANNOVAR | 仅预留表格合并入口 | 软件和 humandb 需手动许可下载；不是 VEP 运行必需项 |
| 新抗原 | 重建突变蛋白并生成任意 8-15mer | 增加表达过滤、转录本优先级、克隆性和质量过滤 |
| HLA typing/binding | HLA*LA G-group分型；MHCflurry/NetMHCpan binding | 用真值样本验证，评估HLA-HD/SpecHLA一致性和WES覆盖限制 |
| CNV | CNVkit支持matched/prebuilt normal；mosdepth仅depth QC | 使用同捕获批次 normal 建 pooled reference，评估 purity/ploidy |
| MSI | MSIsensor-pro paired/tumor-only | 为每套 BED/平台准备位点表和基线，并做阴阳性样本校准 |
| SV | Manta 兼容模块 | Manta 已归档；应评估 GRIDSS2/Delly 等替代及 WES 灵敏度 |
| TMB | VEP consequence + 配对深度/AF + 群体频率严格筛选 | 为panel构建有效编码区并按癌种校准阈值 |
| 报告 | 中文文本 + MultiQC；另有 demo HTML | 建立通用中文 HTML、审计追踪、版本和免责声明 |

## 环境与工具分层

| 环境 | 默认/可选 | 内容 |
|---|---|---|
| `big_wes_pipeline_env` | 默认 | FastQC、fastp、BWA、Samtools、BCFtools、BEDtools、Picard、GATK、SnpEff、mosdepth、MSIsensor-pro、MultiQC |
| `wes_vep_env` | 默认 | Ensembl VEP 115.2 与 Perl/HTSlib 依赖 |
| `wes_hla_env` | 可选 | MHCflurry；模型需要另行下载 |
| `wes_hla_typing_env` | 可选、Linux | HLA*LA 1.0.4；graph需单独下载和prepare |
| `wes_cnv_env` | 可选 | CNVkit 0.9.12，避免 R/Python 依赖影响核心环境 |
| `wes_sv_env` | 可选、Linux | Manta 1.6.0，隔离其旧运行时 |

NetMHCpan、ANNOVAR、COSMIC 不应由自动脚本绕过许可下载。Qualimap 可作为可选
QC 工具，但当前 Samtools/Picard/Mosdepth 已覆盖基础 QC，因此不作为安装阻断项。

## 参考数据分层

### 全局必需

- GRCh38 FASTA、`.fai`、`.dict` 和 BWA index。
- 与 GRCh38 contig 命名一致的捕获 BED。
- VEP 115 GRCh38 cache。
- 新抗原模块使用的 Ensembl 115 `pep.all` 蛋白 FASTA。

### BQSR/体细胞推荐

- dbSNP、Mills indels、1000G known-sites。
- Mutect2 germline population resource。
- Mutect2 common biallelic SNP资源，用于污染估计。
- 同平台、同捕获试剂盒构建的 Panel of Normals。

### 项目或批次特异

- CNVkit pooled normal reference。
- MSIsensor-pro site list 和 tumor-only baseline。
- HLA allele 结果或 HLA 分型数据库。
- 捕获试剂盒有效编码区，用于 TMB 分母校准。

## 已完成的主链加固

1. Mutect2 F1R2自动进入 `LearnReadOrientationModel`，orientation priors传给过滤步骤。
2. common biallelic位点自动执行肿瘤/正常 `GetPileupSummaries` 和 `CalculateContamination`。
3. TMB严格解析VEP CSQ，并输出纳入、排除及原因明细。
4. HLA*LA已作为Linux独立环境接入，保留完整G-group调用并生成binding兼容等位基因。

## 下一轮更新优先级

### P0：让体细胞主链达到可验证状态

1. 增加 PoN、germline resource、common-variants和参考资源的严格preflight及自动索引。
2. 使用SEQC2/内部真值集校准Mutect2过滤、污染估计和低频检出限。
3. 为不同panel建立版本化TMB有效编码区和癌种阈值。

### P1：完善下游结果

1. 把通用中文 HTML 报告接入步骤 13，并展示软件/数据库版本。
2. 使用已知HLA真值样本验证HLA*LA分型，并在报告中展示歧义、覆盖和质量值。
3. CNVkit 使用 matched/pooled normal，并输出 purity/ploidy 相关提示。
4. MSI 为每套捕获 BED 建立版本化位点表和验证阈值。

### P2：队列和工程化

1. 增加 gVCF 与多样本联合分型。
2. 将步骤迁移到 Nextflow/Snakemake，支持调度器、scatter/gather 和容器。
3. 增加小型合成回归测试、真实阳性样本集和 CI。
4. 建立配置 schema、资源校验、版本冻结和结果 provenance。

## 推荐验收标准

- 安装：所有环境创建成功，`ENV_ROOT/manifests` 有请求和解析版本清单。
- 项目生成：单样本和配对模式均能在空目录生成，路径不依赖本机用户名。
- 主流程：根目录一次 `bash run_pipeline.sh` 完成所有已启用模块。
- 恢复：故障后 `from ROLE STEP` 可续跑，已存在 FASTQ/BAM 不重复计算。
- 变异：使用 GIAB/SEQC2 或内部真值样本统计 sensitivity、precision、F1。
- 下游：CNV/MSI/TMB/HLA 各自使用经过版本化和验证的参考资源。
- 报告：明确区分“未运行、测试 fallback、正式结果”，并记录软件和数据库版本。
