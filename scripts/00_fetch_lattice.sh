#!/usr/bin/env bash
# Sparse shallow clone of the FCC-ee LCC_107 Z lattice into data/.
set -euo pipefail
cd "$(dirname "$0")/.."
if [ -f data/acc-models-fcc-ee/lattices/z/fccee_z.madx ]; then
  echo "lattice exists"; exit 0
fi
git clone --depth 1 --branch LCC_107 --filter=blob:none --sparse \
  https://gitlab.cern.ch/acc-models/fcc/acc-models-fcc-ee.git data/acc-models-fcc-ee
git -C data/acc-models-fcc-ee sparse-checkout set lattices/z
