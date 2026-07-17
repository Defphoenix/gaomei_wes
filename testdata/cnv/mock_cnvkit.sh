#!/bin/bash
set -euo pipefail

if [ "${1:-}" = "--version" ] || [ "${1:-}" = "version" ]; then
    echo "CNVkit 0.9.12 mock"
    exit 0
fi

command_name="${1:-}"
shift || true

arg_value() {
    local wanted="$1"
    shift
    while [ $# -gt 1 ]; do
        if [ "$1" = "${wanted}" ]; then
            echo "$2"
            return 0
        fi
        shift
    done
    return 1
}

output=$(arg_value -o "$@")
mkdir -p "$(dirname "${output}")"

case "${command_name}" in
    target)
        input="$1"
        cp "${input}" "${output}"
        ;;
    antitarget)
        printf 'chr1\t1000\t2000\tantitarget\n' > "${output}"
        ;;
    coverage)
        printf 'chromosome\tstart\tend\tgene\tdepth\tlog2\nchr1\t0\t1000\tGENE1\t100\t0\n' > "${output}"
        ;;
    reference)
        printf 'chromosome\tstart\tend\tgene\tlog2\tspread\nchr1\t0\t1000\tGENE1\t0\t0.1\n' > "${output}"
        ;;
    fix)
        printf 'chromosome\tstart\tend\tgene\tlog2\tdepth\tweight\nchr1\t0\t1000\tGENE1\t0.8\t180\t1\n' > "${output}"
        ;;
    segment)
        printf 'chromosome\tstart\tend\tgene\tlog2\tprobes\nchr1\t0\t1000\tGENE1\t0.8\t1\n' > "${output}"
        ;;
    call)
        printf 'chromosome\tstart\tend\tgene\tlog2\tprobes\tcn\nchr1\t0\t1000\tGENE1\t0.8\t1\t3\n' > "${output}"
        ;;
    export)
        export_type="$1"
        if [ "${export_type}" = "seg" ]; then
            printf 'ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\nmock\tchr1\t1\t1000\t1\t0.8\n' > "${output}"
        else
            printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tN\t<DUP>\t.\tPASS\tEND=1000;CN=3\n' > "${output}"
        fi
        ;;
    *)
        echo "Unsupported mock CNVkit command: ${command_name}" >&2
        exit 1
        ;;
esac
