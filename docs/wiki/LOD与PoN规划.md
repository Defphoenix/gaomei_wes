# LOD与PoN规划

v1.0 已完成单个真实配对样本的端到端流程验证。下一阶段应转向队列资源和方法学性能，
不能只以“命令能跑完”作为生产验收标准。

## 1. Panel of Normals

### 样本要求

- 建议起步30份，正式锁定优先50-100份normal。
- 与正式样本使用相同测序平台、捕获试剂盒、读长、参考版本和分析流程。
- 排除肿瘤污染、明显污染、低覆盖和异常文库。
- 尽量覆盖多个批次和建库日期，用于学习稳定的技术伪影。

### 构建原则

```text
normal FASTQ
-> same QC/alignment/dedup/BQSR
-> Mutect2 tumor-only technical calls per normal
-> GenomicsDBImport
-> CreateSomaticPanelOfNormals
-> panel_of_normals.vcf.gz + index
```

PoN必须记录样本列表、流程commit、参考资源版本、capture BED、构建日期和校验值。

## 2. CNV pooled reference

- 使用同平台normal BAM构建CNVkit reference。
- 样本应通过覆盖均一性、性别、污染和批次QC。
- 分别评估单matched normal和pooled reference的噪声、分段稳定性及已知阳性召回。
- capture kit或实验体系变化后重新评估，不默认跨panel共用。

## 3. LOD梯度

建议设计VAF梯度：

```text
20%, 10%, 5%, 2%, 1%, 0.5%
```

每个VAF至少覆盖：

- SNV和短InDel。
- 不同GC、重复区域和比对难度。
- 不同深度，例如100x、200x、500x。
- 至少3个技术重复或独立抽样seed。

统计指标：

| 指标 | 说明 |
|---|---|
| Sensitivity | 预期阳性中被检出的比例 |
| Precision/PPV | 检出结果中真实阳性的比例 |
| FDR | 假阳性比例 |
| VAF bias | 观察VAF相对预期VAF偏差 |
| Repeatability | 批内重复一致性 |
| Reproducibility | 批间/日期/操作者一致性 |

## 4. MSI和TMB

- MSI位点list必须由同一GRCh38和capture BED生成或筛选。
- 使用已知MSS、MSI-L、MSI-H样本建立阈值，不直接照搬其他panel阈值。
- TMB denominator应由有效捕获CDS区域去重合并后计算。
- 使用外部验证集评估不同癌种、纯度和深度下的TMB偏差。

## 5. 变更控制

以下任一变化都应触发再验证或桥接验证：

- 捕获试剂盒版本变化。
- 测序平台、读长或建库方案变化。
- GRCh38参考或contig规范变化。
- GATK/VEP主版本变化。
- PoN、CNV reference、MSI list或TMB BED变化。
- 核心过滤阈值变化。

## 6. 推荐交付物

1. 锁定的软件manifest和Conda explicit lock。
2. 数据库manifest、URL、日期和SHA256。
3. PoN和CNV reference样本清单。
4. LOD、准确性、精密度、特异性和最低肿瘤纯度报告。
5. MSI/TMB阈值报告。
6. 固定的流程tag、配置模板和回归测试结果。

