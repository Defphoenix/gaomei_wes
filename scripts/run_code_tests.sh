#!/bin/bash
# Lightweight regression tests that do not require reference data or conda tools.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
TEST_ROOT="${TMPDIR:-/tmp}/gaomei_wes_code_tests.$$"

cleanup() {
    case "${TEST_ROOT}" in
        "${TMPDIR:-/tmp}"/gaomei_wes_code_tests.*) rm -rf "${TEST_ROOT}" ;;
    esac
}
trap cleanup EXIT
mkdir -p "${TEST_ROOT}"

echo "[1/5] Bash syntax"
bash -n "${PROJECT_DIR}/run_pipeline.sh" "${PROJECT_DIR}/config.sh" "${PROJECT_DIR}"/scripts/*.sh

echo "[2/5] Python syntax"
python3 - "${PROJECT_DIR}/scripts" <<'PY'
import ast
import pathlib
import sys

for path in pathlib.Path(sys.argv[1]).glob("*.py"):
    ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
PY

echo "[3/5] MSI parser stale-output regression"
prefix="${TEST_ROOT}/SH05677_vs_SH05678_msi"
cp "${PROJECT_DIR}/testdata/msi/mock_msisensor_output.tsv" "${prefix}"
printf 'MSI score(%%): 5677\n' > "${prefix}_call_summary.txt"
printf 'sample\tmode\tmsi_status\tmsi_score_percent\nold\tpaired\tMSI-H\t5677\n' > "${prefix}_call.tsv"
python3 "${PROJECT_DIR}/scripts/msi_call_from_msisensor.py" \
    --prefix "${prefix}" \
    --sample SH05677_vs_SH05678 \
    --mode paired \
    --out "${prefix}_call.tsv" \
    --summary "${prefix}_call_summary.txt" \
    --low-threshold 3 \
    --high-threshold 20 >/dev/null
awk -F '\t' 'NR==2 && $3=="MSS" && $4>0.05 && $4<0.06 && $5==3 && $6==5425 {ok=1} END{exit !ok}' "${prefix}_call.tsv"

echo "[4/5] HLA parser"
python3 "${PROJECT_DIR}/scripts/parse_hlala_results.py" \
    --input "${PROJECT_DIR}/testdata/hla/mock_R1_bestguess_G.txt" \
    --output "${TEST_ROOT}/hla.tsv" \
    --alleles-output "${TEST_ROOT}/hla_alleles.txt" \
    --sample mock >/dev/null
test "$(wc -l < "${TEST_ROOT}/hla_alleles.txt" | tr -d ' ')" -ge 1

echo "[5/5] CNVkit matched-normal command regression"
export CNV_TEST_ROOT="${TEST_ROOT}/cnv_case"
bash "${PROJECT_DIR}/run_pipeline.sh" \
    --config "${PROJECT_DIR}/testdata/cnv/config_cnv_matched_test.sh" \
    step 8 >/dev/null
cnv_root="${CNV_TEST_ROOT}/cnv"
test -s "${cnv_root}/mock_tumor_vs_normal.cnr"
test -s "${cnv_root}/mock_tumor_vs_normal.call.cns"
grep -q '^基线类型: single_matched_normal$' "${cnv_root}/mock_tumor_vs_normal_cnv_summary.txt"
grep -q '^  gain: 1$' "${cnv_root}/mock_tumor_vs_normal_cnv_summary.txt"

echo "All lightweight code tests passed."
