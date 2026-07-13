#!/bin/bash
#===============================================================================
# prepare_data.sh - 下载参考基因组、数据库和测试数据
#
# 说明: 本脚本用于下载突变分析流程所需的所有参考文件
#       包括: 参考基因组、BWA索引、已知位点数据库、注释数据库、测试数据
#
# 用法:
#   bash prepare_data.sh              # 下载全部数据
#   bash prepare_data.sh reference    # 仅下载参考基因组
#   bash prepare_data.sh database     # 仅下载数据库
#   bash prepare_data.sh snpeff       # 仅下载SnpEff数据库
#   bash prepare_data.sh vep          # 仅下载VEP缓存
#   bash prepare_data.sh testdata     # 仅下载测试数据
#   bash prepare_data.sh check        # 仅检查已下载文件
#
# 注意: 参考基因组约3.5GB, 全部数据约15-20GB
#       下载时间取决于网络速度
#===============================================================================

set -euo pipefail

#---------------------------------------
# 配置下载参数
#---------------------------------------
# 项目目录
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# 下载目录
REF_DIR="${PIPELINE_DIR}/reference"       # 参考基因组
DB_DIR="${PIPELINE_DIR}/database"         # 数据库
TEST_DIR="${PIPELINE_DIR}/testdata"       # 测试数据
SNPEFF_DIR="${PIPELINE_DIR}/snpeff_db"    # SnpEff数据库
VEP_DIR="${PIPELINE_DIR}/vep_cache"       # VEP缓存

# 下载工具参数 (支持断点续传)
WGET_OPTS="-c -q --show-progress"
CURL_OPTS="-L -C - --progress-bar"

# 是否使用 wget (true) 或 curl (false)
USE_WGET=true

# GATK Resource Bundle 基础URL
GATK_BUNDLE_URL="https://storage.googleapis.com/genomics-publicldata/resources/broad/hg38/v0"

# 1000 Genomes 测试数据URL
G1000_URL="https://storage.googleapis.com/genomics-public-data/1000-genomes"

#---------------------------------------
# 通用下载函数
#---------------------------------------
download_file() {
    local url="$1"
    local output="$2"
    local desc="${3:-文件}"

    if [ -f "${output}" ]; then
        echo "  [已存在] ${desc}: ${output}"
        return 0
    fi

    echo "  [下载中] ${desc}..."
    echo "    URL: ${url}"
    echo "    目标: ${output}"

    mkdir -p "$(dirname "${output}")"

    if [ "${USE_WGET}" = true ]; then
        wget ${WGET_OPTS} -O "${output}" "${url}"
    else
        curl ${CURL_OPTS} -o "${output}" "${url}"
    fi

    if [ $? -eq 0 ] && [ -f "${output}" ]; then
        local size=$(du -h "${output}" | cut -f1)
        echo "  [完成] ${desc} (${size})"
    else
        echo "  [失败] ${desc} 下载失败!"
        return 1
    fi
}

#---------------------------------------
# 打印信息
#---------------------------------------
print_header() {
    echo ""
    echo "════════════════════════════════════════════════════════════"
    echo "  $1"
    echo "════════════════════════════════════════════════════════════"
    echo ""
}

print_step() {
    echo "────────────────────────────────────────────────────────"
    echo "  >> $1"
    echo "────────────────────────────────────────────────────────"
}

#=======================================
#  第一部分: 参考基因组
#=======================================
download_reference() {
    print_header "下载参考基因组 (hg38)"

    mkdir -p "${REF_DIR}"

    #---------------------------------------
    # 1.1 参考基因组 FASTA
    #---------------------------------------
    print_step "1.1 参考基因组 FASTA (~900MB)"
    download_file \
        "${GATK_BUNDLE_URL}/Homo_sapiens_assembly38.fasta.gz" \
        "${REF_DIR}/hg38.fa.gz" \
        "参考基因组 (hg38 FASTA)"

    # 解压
    if [ -f "${REF_DIR}/hg38.fa.gz" ] && [ ! -f "${REF_DIR}/hg38.fa" ]; then
        echo "  解压参考基因组..."
        gunzip -k "${REF_DIR}/hg38.fa.gz"
    fi

    #---------------------------------------
    # 1.2 参考基因组索引 (.fai)
    #---------------------------------------
    print_step "1.2 参考基因组索引 (.fai)"
    download_file \
        "${GATK_BUNDLE_URL}/Homo_sapiens_assembly38.fasta.fai" \
        "${REF_DIR}/hg38.fa.fai" \
        "参考基因组索引 (.fai)"

    #---------------------------------------
    # 1.3 参考基因组 dict (.dict)
    #---------------------------------------
    print_step "1.3 参考基因组 dict (.dict)"
    download_file \
        "${GATK_BUNDLE_URL}/Homo_sapiens_assembly38.dict" \
        "${REF_DIR}/hg38.dict" \
        "参考基因组 dict (.dict)"

    #---------------------------------------
    # 1.4 BWA索引
    #---------------------------------------
    print_step "1.4 BWA索引 (~3.5GB, 需要较长时间)"
    echo "  使用bwa index构建索引 (如已有则跳过)..."

    if [ ! -f "${REF_DIR}/hg38.fa.bwt" ]; then
        if command -v bwa &>/dev/null; then
            echo "  运行: bwa index ${REF_DIR}/hg38.fa"
            bwa index "${REF_DIR}/hg38.fa"
            echo "  BWA索引构建完成"
        else
            echo "  [跳过] bwa未安装，请手动运行: bwa index ${REF_DIR}/hg38.fa"
            echo "  或从以下URL下载预构建索引:"
            echo "    ${GATK_BUNDLE_URL}/Homo_sapiens_assembly38.fasta.bwt"
            echo "    ${GATK_BUNDLE_URL}/Homo_sapiens_assembly38.fasta.pac"
            echo "    ${GATK_BUNDLE_URL}/Homo_sapiens_assembly38.fasta.ann"
            echo "    ${GATK_BUNDLE_URL}/Homo_sapiens_assembly38.fasta.amb"
            echo "    ${GATK_BUNDLE_URL}/Homo_sapiens_assembly38.fasta.sa"
            # 尝试下载预构建索引
            for ext in bwt pac ann amb sa; do
                download_file \
                    "${GATK_BUNDLE_URL}/Homo_sapiens_assembly38.fasta.${ext}" \
                    "${REF_DIR}/hg38.fa.${ext}" \
                    "BWA索引 .${ext}" || true
            done
        fi
    else
        echo "  [已存在] BWA索引已完成"
    fi

    echo ""
    echo "  参考基因组下载完成!"
    echo "  目录: ${REF_DIR}"
    ls -lh "${REF_DIR}"/hg38.fa* 2>/dev/null | awk '{print "    "$NF" ("$5")"}'
}

#=======================================
#  第二部分: 数据库
#=======================================
download_databases() {
    print_header "下载已知位点数据库"

    mkdir -p "${DB_DIR}"

    #---------------------------------------
    # 2.1 dbSNP (必需 - BQSR和注释用)
    #---------------------------------------
    print_step "2.1 dbSNP数据库 (~500MB)"
    download_file \
        "${GATK_BUNDLE_URL}/dbsnp_146.hg38.vcf.gz" \
        "${DB_DIR}/dbsnp_146.hg38.vcf.gz" \
        "dbSNP v146"

    download_file \
        "${GATK_BUNDLE_URL}/dbsnp_146.hg38.vcf.gz.tbi" \
        "${DB_DIR}/dbsnp_146.hg38.vcf.gz.tbi" \
        "dbSNP索引"

    #---------------------------------------
    # 2.2 已知Indels (BQSR用)
    #---------------------------------------
    print_step "2.2 已知Indels数据库"
    download_file \
        "${GATK_BUNDLE_URL}/Homo_sapiens_assembly38.known_indels.vcf.gz" \
        "${DB_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz" \
        "已知Indels"

    download_file \
        "${GATK_BUNDLE_URL}/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi" \
        "${DB_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi" \
        "已知Indels索引"

    #---------------------------------------
    # 2.3 HapMap (VQSR用, 可选)
    #---------------------------------------
    print_step "2.3 HapMap数据库"
    download_file \
        "${GATK_BUNDLE_URL}/hapmap_3.3.hg38.vcf.gz" \
        "${DB_DIR}/hapmap_3.3.hg38.vcf.gz" \
        "HapMap 3.3"

    download_file \
        "${GATK_BUNDLE_URL}/hapmap_3.3.hg38.vcf.gz.tbi" \
        "${DB_DIR}/hapmap_3.3.hg38.vcf.gz.tbi" \
        "HapMap索引"

    #---------------------------------------
    # 2.4 1000 Genomes Omni (VQSR用, 可选)
    #---------------------------------------
    print_step "2.4 1000G Omni数据库"
    download_file \
        "${GATK_BUNDLE_URL}/1000G_omni2.5.hg38.vcf.gz" \
        "${DB_DIR}/1000G_omni2.5.hg38.vcf.gz" \
        "1000G Omni 2.5"

    download_file \
        "${GATK_BUNDLE_URL}/1000G_omni2.5.hg38.vcf.gz.tbi" \
        "${DB_DIR}/1000G_omni2.5.hg38.vcf.gz.tbi" \
        "1000G Omni索引"

    #---------------------------------------
    # 2.5 1000 Genomes Phase1 SNPs (VQSR用, 可选)
    #---------------------------------------
    print_step "2.5 1000G Phase1 SNPs"
    download_file \
        "${GATK_BUNDLE_URL}/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
        "${DB_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz" \
        "1000G Phase1 SNPs"

    download_file \
        "${GATK_BUNDLE_URL}/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi" \
        "${DB_DIR}/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi" \
        "1000G Phase1索引"

    #---------------------------------------
    # 2.6 Mills and 1000G Indels (BQSR用)
    #---------------------------------------
    print_step "2.6 Mills & 1000G Gold Standard Indels"
    download_file \
        "${GATK_BUNDLE_URL}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
        "${DB_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
        "Mills & 1000G Indels"

    download_file \
        "${GATK_BUNDLE_URL}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" \
        "${DB_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi" \
        "Mills Indels索引"

    #---------------------------------------
    # 2.7 外显子组目标区域BED (示例 - Agilent V6)
    #---------------------------------------
    print_step "2.7 外显子组目标区域BED (示例)"
    echo "  注意: 请根据实际使用的捕获试剂盒替换此文件!"
    echo "  这里使用一个通用的编码区BED作为示例"

    # 从UCSC下载hg38编码区作为示例
    if [ ! -f "${DB_DIR}/target_regions.bed" ]; then
        echo "  生成示例目标区域BED (chr21: 最小染色体用于测试)..."
        # 创建一个简单的chr21测试区域BED
        cat > "${DB_DIR}/target_regions.bed" << 'BED_EOF'
chr21	9000000	9500000	test_region_1
chr21	10000000	10500000	test_region_2
chr21	11000000	11500000	test_region_3
chr21	12000000	12500000	test_region_4
chr21	13000000	13500000	test_region_5
chr21	14000000	14500000	test_region_6
chr21	15000000	15500000	test_region_7
chr21	16000000	16500000	test_region_8
BED_EOF
        echo "  示例BED已生成: ${DB_DIR}/target_regions.bed"
    fi

    echo ""
    echo "  数据库下载完成!"
    echo "  目录: ${DB_DIR}"
    ls -lh "${DB_DIR}"/*.vcf.gz 2>/dev/null | awk '{print "    "$NF" ("$5")"}'
}

#=======================================
#  第三部分: SnpEff 数据库
#=======================================
download_snpeff_db() {
    print_header "下载SnpEff注释数据库"

    mkdir -p "${SNPEFF_DIR}"

    if command -v snpEff &>/dev/null; then
        echo "  使用snpEff download命令下载..."
        snpEff download -dataDir "${SNPEFF_DIR}" GRCh38.105
    else
        echo "  SnpEff未安装，尝试直接下载数据库..."
        echo ""
        echo "  方法1: 安装SnpEff后运行:"
        echo "    snpEff download -dataDir ${SNPEFF_DIR} GRCh38.105"
        echo ""
        echo "  方法2: 手动下载:"
        echo "    URL: https://snpeff.blob.core.windows.net/databases/v5_2/GRCh38.105/snpEff_databases.tgz"
        echo "    解压到: ${SNPEFF_DIR}"
        echo ""

        # 尝试直接下载
        download_file \
            "https://snpeff.blob.core.windows.net/databases/v5_2/GRCh38.105/snpEff_databases.tgz" \
            "${SNPEFF_DIR}/snpEff_databases.tgz" \
            "SnpEff GRCh38.105数据库" || \
        echo "  自动下载失败，请手动下载"
    fi
}

#=======================================
#  第四部分: VEP 缓存
#=======================================
download_vep_cache() {
    print_header "下载VEP注释缓存"

    mkdir -p "${VEP_DIR}"

    if command -v vep &>/dev/null; then
        echo "  使用VEP安装器下载缓存..."
        vep_install \
            --CACHEDIR "${VEP_DIR}" \
            --SPECIES homo_sapiens \
            --ASSEMBLY GRCh38 \
            --CACHE_VERSION 110 \
            --AUTO cf
    else
        echo "  VEP未安装，请手动安装后下载:"
        echo ""
        echo "  方法1: conda安装"
        echo "    conda install -c bioconda ensembl-vep"
        echo "    vep_install --CACHEDIR ${VEP_DIR} --AUTO cf"
        echo ""
        echo "  方法2: 手动下载"
        echo "    URL: https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz"
        echo "    解压到: ${VEP_DIR}"
        echo ""

        download_file \
            "https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz" \
            "${VEP_DIR}/homo_sapiens_vep_110_GRCh38.tar.gz" \
            "VEP GRCh38缓存" || \
        echo "  自动下载失败，请手动下载"
    fi
}

#=======================================
#  第五部分: 测试数据 (已修改为直接从本地提取和模拟)
#=======================================
download_testdata() {
    print_header "生成测试数据 (从本地hg38提取chr21并模拟reads)"

    mkdir -p "${TEST_DIR}"

    # 指定你真实的本地参考基因组绝对路径
    local REAL_REF="/PUBLIC/gomics/guofenghua/project/project_01_wespipline/pipline/reference_data/hg38/Homo_sapiens_assembly38.fasta"
    local REGION="chr21:9000000-12000000"

    echo "  测试策略: 使用本地参考基因组截取 ${REGION} 区域"
    echo "  并直接使用 wgsim 生成双端模拟测序数据"
    echo ""

    # 检查本地参考基因组是否存在
    if [ ! -f "${REAL_REF}" ]; then
        echo "  [错误] 找不到本地参考基因组: ${REAL_REF}"
        echo "  请检查路径是否正确！"
        return 1
    fi

    #---------------------------------------
    # 1. 提取参考基因组子集
    #---------------------------------------
    print_step "1. 提取 ${REGION} 作为测试参考序列"
    samtools faidx "${REAL_REF}" "${REGION}" > "${TEST_DIR}/test_ref.fa"
    
    if [ ! -s "${TEST_DIR}/test_ref.fa" ]; then
        echo "  [错误] 提取子集失败！请检查 samtools 是否可用。"
        return 1
    fi
    echo "  提取完成: ${TEST_DIR}/test_ref.fa"

    #---------------------------------------
    # 2. 使用 wgsim 模拟 Reads 并压缩
    #---------------------------------------
    print_step "2. 使用 wgsim 模拟 10,000 条双端测序 Reads"
    if command -v wgsim &>/dev/null; then
        wgsim -e 0.001 -r 0.01 -1 150 -2 150 -d 400 -s 50 -N 10000 \
            "${TEST_DIR}/test_ref.fa" \
            "${TEST_DIR}/test_R1.fastq" \
            "${TEST_DIR}/test_R2.fastq"

        echo "  使用 pigz (多线程) 压缩 FASTQ 文件..."
        pigz -f "${TEST_DIR}/test_R1.fastq"
        pigz -f "${TEST_DIR}/test_R2.fastq"

        echo "  生成成功: test_R1.fastq.gz, test_R2.fastq.gz"
    else
        echo "  [错误] 未找到 wgsim 命令，请确保已在 Conda 环境中！"
        return 1
    fi

    #---------------------------------------
    # 3. 生成测试用目标区域 BED 文件
    #---------------------------------------
    print_step "3. 生成配套的测试目标区域 BED (test_target.bed)"
    cat > "${TEST_DIR}/test_target.bed" << 'BED_EOF'
chr21	9000000	9200000	target_1
chr21	9300000	9500000	target_2
chr21	9600000	9800000	target_3
chr21	10000000	10200000	target_4
chr21	10300000	10500000	target_5
chr21	10600000	10800000	target_6
chr21	11000000	11200000	target_7
chr21	11300000	11500000	target_8
chr21	11600000	11800000	target_9
chr21	12000000	12200000	target_10
BED_EOF
    echo "  测试 BED 已生成: ${TEST_DIR}/test_target.bed"

    echo ""
    echo "  ========================================"
    echo "  所有本地测试数据已成功生成！可以直接开跑流程！"
    echo "  目录: ${TEST_DIR}"
    ls -lh "${TEST_DIR}"/test_* 2>/dev/null | awk '{print "    "$NF" ("$5")"}'
}

#=======================================
#  检查已下载文件
#=======================================
check_files() {
    print_header "检查已下载文件"

    local all_ok=true

    # 参考基因组
    print_step "参考基因组"
    for f in "hg38.fa" "hg38.fa.fai" "hg38.dict"; do
        if [ -f "${REF_DIR}/${f}" ]; then
            echo "  [OK] ${f}"
        else
            echo "  [缺失] ${f}"
            all_ok=false
        fi
    done
    # BWA索引
    for ext in bwt pac ann amb sa; do
        if [ -f "${REF_DIR}/hg38.fa.${ext}" ]; then
            echo "  [OK] hg38.fa.${ext}"
        else
            echo "  [缺失] hg38.fa.${ext}"
            all_ok=false
        fi
    done

    # 数据库
    print_step "数据库"
    for f in "dbsnp_146.hg38.vcf.gz" "Homo_sapiens_assembly38.known_indels.vcf.gz" \
             "hapmap_3.3.hg38.vcf.gz" "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"; do
        if [ -f "${DB_DIR}/${f}" ]; then
            echo "  [OK] ${f}"
        else
            echo "  [缺失] ${f}"
            all_ok=false
        fi
    done

    # SnpEff
    print_step "SnpEff数据库"
    if [ -d "${SNPEFF_DIR}/GRCh38.105" ] || [ -f "${SNPEFF_DIR}/snpEff_databases.tgz" ]; then
        echo "  [OK] SnpEff GRCh38.105"
    else
        echo "  [缺失] SnpEff GRCh38.105"
        all_ok=false
    fi

    # VEP
    print_step "VEP缓存"
    if [ -d "${VEP_DIR}/homo_sapiens" ]; then
        echo "  [OK] VEP GRCh38缓存"
    else
        echo "  [缺失] VEP GRCh38缓存"
        all_ok=false
    fi

    # 测试数据
    print_step "测试数据"
    if ls "${TEST_DIR}"/*fastq* &>/dev/null || ls "${TEST_DIR}"/*bam* &>/dev/null; then
        echo "  [OK] 测试数据已准备"
    else
        echo "  [缺失] 测试数据未准备"
        all_ok=false
    fi

    echo ""
    if [ "${all_ok}" = true ]; then
        echo "  所有文件检查通过!"
    else
        echo "  部分文件缺失，请运行对应的下载命令"
    fi
}

#=======================================
#  生成config.sh更新提示
#=======================================
print_config_tips() {
    print_header "配置提示: 修改config.sh"

    echo "  下载完成后，请修改 config.sh 中的以下路径:"
    echo ""
    echo "  # 参考基因组"
    echo "  REFERENCE_GENOME=\"${REF_DIR}/hg38.fa\""
    echo "  REFERENCE_DICT=\"${REF_DIR}/hg38.dict\""
    echo "  REFERENCE_BWA_INDEX=\"${REF_DIR}/hg38.fa\""
    echo ""
    echo "  # 数据库"
    echo "  DBSNP_VCF=\"${DB_DIR}/dbsnp_146.hg38.vcf.gz\""
    echo "  KNOWN_INDELS_VCF=\"${DB_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz\""
    echo "  MILLS_VCF=\"${DB_DIR}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz\""
    echo ""
    echo "  # SnpEff"
    echo "  SNPEFF_DATA_DIR=\"${SNPEFF_DIR}\""
    echo "  SNPEFF_CONFIG=\"${SNPEFF_DIR}/snpEff.config\""
    echo ""
    echo "  # VEP"
    echo "  VEP_CACHE_DIR=\"${VEP_DIR}\""
    echo ""
    echo "  # 目标区域"
    echo "  INTERVAL_FILE=\"${TEST_DIR}/test_target.bed\""
    echo ""
    echo "  # 测试数据"
    echo "  FASTQ_R1=\"${TEST_DIR}/test_R1.fastq.gz\""
    echo "  FASTQ_R2=\"${TEST_DIR}/test_R2.fastq.gz\""
    echo ""
}

#=======================================
#  主入口
#=======================================
main() {
    local action="${1:-all}"

    case "${action}" in
        reference|ref)
            download_reference
            ;;
        database|db)
            download_databases
            ;;
        snpeff)
            download_snpeff_db
            ;;
        vep)
            download_vep_cache
            ;;
        testdata|test)
            download_testdata
            ;;
        check)
            check_files
            ;;
        all|full)
            download_reference
            download_databases
            download_snpeff_db
            download_vep_cache
            download_testdata
            check_files
            print_config_tips
            ;;
        help|--help|-h)
            echo "用法: bash $0 [选项]"
            echo ""
            echo "选项:"
            echo "  all (默认)    下载全部数据"
            echo "  reference     仅参考基因组"
            echo "  database      仅数据库"
            echo "  snpeff        仅SnpEff数据库"
            echo "  vep           仅VEP缓存"
            echo "  testdata      仅测试数据"
            echo "  check         检查已下载文件"
            ;;
        *)
            echo "未知操作: ${action}"
            echo "运行 bash $0 help 查看帮助"
            exit 1
            ;;
    esac
}

main "$@"
