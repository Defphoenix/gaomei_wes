# Changelog

## v1.0.0 - 2026-07-17

首个完整研发版，支持单样本胚系和 tumor-normal 配对 WES。主要能力：

- 一键安装的 prefix-based mamba 环境及精确软件 manifest。
- FASTQ QC、修剪、BWA、去重、比对后 QC、BQSR 和 BAM 完整性检查。
- 配对项目中 normal/tumor 仅运行到 BQSR/HLA，Mutect2 只在 somatic 阶段运行一次。
- Mutect2 F1R2 方向偏倚模型、可选污染估计、PoN 和 AF-only germline resource。
- VEP 115 离线注释、8-15mer 新抗原、HLA*LA 分型和 MHCflurry/NetMHCpan binding。
- CNVkit matched/prebuilt normal CNV；mosdepth 无基线时仅输出 depth QC。
- MSIsensor-pro paired MSI、严格 VEP TMB、覆盖度和中文汇总报告。
- 单步、断点、状态检查和真实配对样本回归测试。

真实样本验证期间修复：

- GATK Java 17、SnpEff Java 21 环境隔离。
- Picard 3.4 `MINIMUM_PCT` 参数兼容。
- `samtools markdup` 参数兼容。
- BQSR 输出未完成时被下游提前读取的问题。
- VEP cache 生效配置和空 VCF 处理。
- MSI 解析器读取自身旧摘要导致错误 MSI-H 的问题。
- Mutect2 测试下采样值过低的问题。
- mosdepth depth-ratio 被误标为正式 CNV 的问题。
- `bedtools multicov` 在大型 WES BED 上重复且耗时的问题。

v1.0.0 仍是研发/验证版。正式临床或工业发布前需完成 PoN、gnomAD AF-only、
同平台 CNV reference、panel MSI/TMB 阈值、LOD/重复性/准确性和临床知识库验证。
