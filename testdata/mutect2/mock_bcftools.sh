#!/bin/bash
set -euo pipefail

case "$1" in
    view)
        output=""
        input="${!#}"
        shift
        while [ $# -gt 1 ]; do
            if [ "$1" = "-o" ]; then
                output="$2"
                break
            fi
            shift
        done
        cp "${input}" "${output}"
        ;;
    index)
        input="${!#}"
        touch "${input}.tbi"
        ;;
    *)
        echo "Unsupported mock bcftools invocation: $*" >&2
        exit 1
        ;;
esac

