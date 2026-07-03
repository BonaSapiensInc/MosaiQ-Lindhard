#!/usr/bin/env bash
# Build a minimal APS submission zip: two-fermi.tex + cited figure PDFs only.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MANUSCRIPT="$ROOT/manuscript"
FIGURES="$MANUSCRIPT/figures"
OUTPUT="$MANUSCRIPT/submission.zip"
TEX="$MANUSCRIPT/two-fermi.tex"

if [[ ! -f "$TEX" ]]; then
  echo "error: missing $TEX" >&2
  exit 1
fi

# Figure stems cited in two-fermi.tex (keep in sync with \panelpagefigure / \includegraphics).
CITED=(
  chi_e_re chi_e_im chi_i_re chi_i_im chi_ei_re chi_ei_im
  S_ee_contour S_ii_contour S_ei_contour
  plasmon_dispersion
  t0_analytic_re_contour t0_analytic_im_contour t0_limit_validation
  thermal_anatomy_a thermal_anatomy_b thermal_anatomy_c thermal_anatomy_d
  Sq_gamma_ee Sq_gamma_ii Sq_gamma_ei
)

missing=()
for stem in "${CITED[@]}"; do
  if [[ ! -f "$FIGURES/${stem}.pdf" ]]; then
    missing+=("${stem}.pdf")
  fi
done
if ((${#missing[@]} > 0)); then
  echo "error: missing cited figure(s): ${missing[*]}" >&2
  exit 1
fi

find "$MANUSCRIPT" -name '.DS_Store' -delete 2>/dev/null || true

staging="$(mktemp -d)"
trap 'rm -rf "$staging"' EXIT
mkdir -p "$staging/figures"
cp "$TEX" "$staging/"
for stem in "${CITED[@]}"; do
  cp "$FIGURES/${stem}.pdf" "$staging/figures/"
done

rm -f "$OUTPUT"
(
  cd "$staging"
  zip -r "$OUTPUT" two-fermi.tex figures/*.pdf
)

echo "Created $OUTPUT"
echo "Contents:"
unzip -l "$OUTPUT"
