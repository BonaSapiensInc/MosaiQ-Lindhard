#!/usr/bin/env bash
# Regenerate static S(q) data and manuscript figures (run from repository root).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

SIM="./simulator/build/mosaiq_simulator"
if [[ ! -x "$SIM" ]]; then
  echo "Error: build the simulator first:" >&2
  echo "  cmake -S simulator -B simulator/build -DCMAKE_BUILD_TYPE=Release" >&2
  echo "  cmake --build simulator/build --target mosaiq_simulator" >&2
  exit 1
fi

RS="${1:-1}"
GAMMAS="${2:-10,50,100,150}"

echo "Running gamma sweep: r_s=${RS}, Gamma=${GAMMAS} (q in [0.1, 25.0], step 0.05)..."
echo "Note: large q and dynamic omega bounds make this slow (often several hours)."
"$SIM" --gamma-sweep "$RS" "$GAMMAS"

echo "Rendering Sq_gamma PDF figures..."
python3 scripts/plot_Sq_gamma.py

echo "Done. Figures: output/Sq_gamma_{ee,ii,ei}.pdf"
