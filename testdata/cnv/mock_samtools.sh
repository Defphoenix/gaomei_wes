#!/bin/bash
set -euo pipefail

case "${1:-}" in
    --version) echo "samtools mock" ;;
    quickcheck) exit 0 ;;
    index) : > "${2}.bai" ;;
    *) exit 0 ;;
esac
