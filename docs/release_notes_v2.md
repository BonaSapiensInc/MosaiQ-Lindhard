# MosaiQ-Lindhard v2.0.0-dual-pathway

**Release tag:** `v2.0.0-dual-pathway`  
**Previous:** `v1.0.2-stratonovich-update`  
**Archive:** [Zenodo DOI 10.5281/zenodo.21331277](https://doi.org/10.5281/zenodo.21331277)

This major release freezes the dual-pathway strong-coupling architecture in the public engine and documentation, alongside the deterministic complex polylogarithm core used by Pathway~A.

---

## Major feature — Dual-pathway strong-coupling architecture

Two parallel frameworks share the same causal Lindhard input $\chi^L$ and thermodynamic distance $f(\Gamma,r_s,\tau)$, but play distinct roles (manuscript Sec.~V / Appendix):

| Pathway | Role | Engine status |
|---------|------|---------------|
| **A — PolyLog-RPA** | Formal scalar diagrammatic definition $\chi^{(s)}=\chi^L\mathrm{Li}_s(v\chi^L)$ with $s=f$ | Live scalar diagnostic (`--pathway polylog-rpa --scalar-diagnostic`) |
| **B — LFC–Zeta / Zeta-RPA** | Production local-field-correction skeleton $\chi^\zeta=\chi^L/(1-W_\zeta v\chi^L)$ | Live multi-component default (`--pathway zeta-rpa`) |

Pathway~B remains the production route for all manuscript multi-component numerics. Scalar PolyLog does **not** promote to the coupled two-component RPA tensor in this release.

---

## Core math — Complex polylogarithm $\mathrm{Li}_s(z)$

New deterministic evaluator (`ComplexPolylogarithm`):

- **Fast path:** power series for $|z| < 0.75$
- **Heavy path:** Mellin integral on the real line after $t=e^u$, summed by truncated **sinc / trapezoidal** quadrature (analytic continuation for $|z|\ge 0.75$)
- **Boundary:** $\mathrm{Li}_s(1)=\zeta(s)$ via Borwein for $s>1$; real branch cut $z\ge 1$ refused except the exact unit point
- Domain extended to **$s>0$** so Pathway~A can set $s=f$ through weak and intermediate coupling

---

## Engine upgrades

- **`PolyLogRPA` module:** scalar Pathway~A evaluator with weak-coupling identity gate ($s\to 0$ → ordinary RPA)
- Shared coupling shape factor $f$ with Zeta-RPA weights
- CLI pathway token `polylog-rpa` (alias `polylog`); multi-component / $\Gamma$-sweep modes remain on Pathway~B
- Plasmon pole extractor routes dielectric roots for Standard RPA, PolyLog-RPA, and Zeta-RPA

---

## Validation — 3-way plasmon dispersion cross-check

`--scalar-diagnostic` now extracts poles for **Standard RPA**, **PolyLog-RPA**, and **Zeta-RPA** at each $\bar{q}$ into `output/output_zeta_rpa_dispersion.dat`.

Visualization:

```bash
./simulator/build/mosaiq_simulator --scalar-diagnostic --gamma 100 2.0 10000
python3 scripts/plot_pathway_ab_dispersion.py
# → output/pathway_ab_comparison.pdf
```

---

## See also

- Theory: [`docs/theory.md`](theory.md)
- CLI: [`docs/usage.md`](usage.md)
- Zeta integration blueprint: [`docs/zeta_rpa_integration.md`](zeta_rpa_integration.md)
- Manuscript: [`manuscript/two-fermi.tex`](../manuscript/two-fermi.tex)
