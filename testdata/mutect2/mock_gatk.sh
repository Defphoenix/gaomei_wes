#!/bin/bash
set -euo pipefail

tool=""
for value in "$@"; do
    case "${value}" in
        LearnReadOrientationModel|GetPileupSummaries|CalculateContamination|FilterMutectCalls)
            tool="${value}"
            break
            ;;
    esac
done

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

case "${tool}" in
    LearnReadOrientationModel)
        output=$(arg_value -O "$@")
        printf 'mock orientation model\n' > "${output}"
        ;;
    GetPileupSummaries)
        output=$(arg_value -O "$@")
        printf 'contig\tposition\tref_count\talt_count\nchr1\t120\t80\t20\n' > "${output}"
        ;;
    CalculateContamination)
        output=$(arg_value -O "$@")
        segments=$(arg_value --tumor-segmentation "$@")
        printf 'sample\tcontamination\terror\nmock_pair\t0.012\t0.001\n' > "${output}"
        printf 'contig\tstart\tend\tminor_allele_fraction\nchr1\t1\t1000\t0.49\n' > "${segments}"
        ;;
    FilterMutectCalls)
        output=$(arg_value -O "$@")
        printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t120\t.\tA\tT\t100\tPASS\t.\n' > "${output}"
        printf '%s\n' "$*" > "${output}.command.txt"
        ;;
    *)
        echo "Unsupported mock GATK invocation: $*" >&2
        exit 1
        ;;
esac

