#!/bin/bash
# Extract, link and index the separately distributed HLA*LA graph package.

set -euo pipefail

ARCHIVE=""
REFERENCE_DIR=""
ENV_PREFIX=""
GRAPH_NAME="PRG_MHC_GRCh38_withIMGT"
EXPECTED_MD5="525a8aa0c7f357bf29fe2c75ef1d477d"

usage() {
    cat <<EOF
Usage: bash scripts/prepare_hlala_graph.sh --archive FILE --reference-dir DIR --env-prefix DIR

  --archive FILE       Downloaded PRG_MHC_GRCh38_withIMGT.tar.gz.
  --reference-dir DIR  WES reference root; graph is stored below DIR/hla/.
  --env-prefix DIR     wes_hla_typing_env prefix.
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --archive) ARCHIVE="$2"; shift 2 ;;
        --reference-dir) REFERENCE_DIR="$2"; shift 2 ;;
        --env-prefix) ENV_PREFIX="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

[ -f "${ARCHIVE}" ] || { echo "Graph archive not found: ${ARCHIVE}" >&2; exit 1; }
[ -d "${REFERENCE_DIR}" ] || { echo "Reference directory not found: ${REFERENCE_DIR}" >&2; exit 1; }
[ -f "${ENV_PREFIX}/conda-meta/history" ] || { echo "HLA typing env not found: ${ENV_PREFIX}" >&2; exit 1; }

if command -v md5sum >/dev/null 2>&1; then
    actual_md5=$(md5sum "${ARCHIVE}" | awk '{print $1}')
elif command -v md5 >/dev/null 2>&1; then
    actual_md5=$(md5 -q "${ARCHIVE}")
else
    actual_md5=""
    echo "[WARN] md5 tool unavailable; archive checksum was not verified" >&2
fi
if [ -n "${actual_md5}" ] && [ "${actual_md5}" != "${EXPECTED_MD5}" ]; then
    echo "Graph archive MD5 mismatch: expected ${EXPECTED_MD5}, got ${actual_md5}" >&2
    exit 1
fi

graph_parent="${REFERENCE_DIR}/hla"
graph_dir="${graph_parent}/${GRAPH_NAME}"
mkdir -p "${graph_parent}"
if [ ! -d "${graph_dir}" ]; then
    echo "[EXTRACT] ${ARCHIVE} -> ${graph_parent}"
    tar -xzf "${ARCHIVE}" -C "${graph_parent}"
fi
[ -d "${graph_dir}" ] || { echo "Archive did not create ${graph_dir}" >&2; exit 1; }

hla_root="${ENV_PREFIX}/opt/hla-la"
hla_graph_link="${hla_root}/graphs/${GRAPH_NAME}"
hla_prepare_bin="${hla_root}/bin/HLA-LA"
mkdir -p "${hla_root}/graphs"
if [ ! -e "${hla_graph_link}" ]; then
    ln -s "${graph_dir}" "${hla_graph_link}"
elif [ "$(cd "${hla_graph_link}" && pwd)" != "$(cd "${graph_dir}" && pwd)" ]; then
    echo "Existing HLA*LA graph path points elsewhere: ${hla_graph_link}" >&2
    exit 1
fi
[ -x "${hla_prepare_bin}" ] || { echo "HLA*LA prepare binary not found: ${hla_prepare_bin}" >&2; exit 1; }

if [ -f "${graph_dir}/serializedGRAPH" ]; then
    echo "[SKIP] graph is already prepared: ${graph_dir}"
else
    echo "[PREPARE] This can take hours and require approximately 40 GB RAM."
    "${hla_prepare_bin}" --action prepareGraph --PRG_graph_dir "${graph_dir}"
fi

echo "HLA*LA graph ready: ${graph_dir}"
echo "HLA*LA graph name: ${GRAPH_NAME}"

