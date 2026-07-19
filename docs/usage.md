# MosaiQ-Lindhard Simulator — Usage

**Binary (typical):** `simulator/build/mosaiq_simulator`  
**Theory:** [`theory.md`](theory.md) (Pathway~A = PolyLog-RPA, Pathway~B = LFC–Zeta)  
**Zeta engine blueprint:** [`zeta_rpa_integration.md`](zeta_rpa_integration.md)

---

## Response pathways (`--pathway`)

| Token | Manuscript role | Engine status |
|-------|-----------------|---------------|
| `zeta-rpa` | **Pathway~B** — LFC–Zeta production (multi-component) | **Live default** |
| `standard-rpa` | Undressed two-component RPA | Live |
| `zeta-rpa-experimental` | Alternate locked weight (diagnostics) | Live |
| `polylog-rpa` | **Pathway~A** — PolyLog-RPA $\chi^L\mathrm{Li}_s(v\chi^L)$, $s=f$ | **Live** (scalar only; requires `--scalar-diagnostic`) |

Aliases: `rpa` → `standard-rpa`, `polylog` → `polylog-rpa`.

```text
--pathway standard-rpa | zeta-rpa | zeta-rpa-experimental | polylog-rpa
            (default: zeta-rpa — Pathway B, multi-component manuscript pipeline)
```

Zeta-RPA is **not** removed when PolyLog-RPA lands. Both pathways remain selectable for side-by-side comparison.

---

## Examples

```bash
# Pathway B — standard LFC–Zeta production (current default)
./simulator/build/mosaiq_simulator 1.0 10000
./simulator/build/mosaiq_simulator --pathway zeta-rpa 1.0 10000

# Legacy undressed RPA
./simulator/build/mosaiq_simulator --pathway standard-rpa 1.0 10000

# Pathway A — PolyLog-RPA (scalar diagnostic; multi-component S(q,ω) stays on zeta-rpa)
./simulator/build/mosaiq_simulator --pathway polylog-rpa --scalar-diagnostic --gamma 50 1.0 10000

# Scalar Zeta diagnostic grids (Pathway B)
./simulator/build/mosaiq_simulator --scalar-diagnostic --pathway zeta-rpa --gamma 50 1.0 1000
# → output/output_zeta_rpa_scalar.dat
# → output/output_zeta_rpa_dispersion.dat

# Static S(q) Γ-sweep (inherits --pathway; default zeta-rpa)
./simulator/build/mosaiq_simulator --gamma-sweep 1 10,50,100,150
./simulator/build/mosaiq_simulator --pathway standard-rpa --gamma-sweep 1 10,50,100,150
```

Stderr prints a pathway header, e.g.  
`ResponsePathway: ZetaRPA (production multi-component Zeta-RPA; strong-coupling window)`.

---

## Performance notes

| Pathway | Cost driver | Recommendation |
|---------|-------------|----------------|
| **B — LFC–Zeta** | One scalar $W_\zeta(f)$ (or channel-wise $W_\zeta^{st}$) then RPA algebra | Production meshes, manuscript figures, Γ-sweeps |
| **A — PolyLog-RPA** | Complex $\mathrm{Li}_s(z)$ at every causal $(q,\omega)$ node | Diagnostic benchmarks and theoretical validation; expect substantially higher wall time |

PolyLog-RPA evaluates complex $\mathrm{Li}_s(z)$ (series for $|z|<0.75$, sinc–Mellin otherwise) at every causal $(q,\omega)$ node. Prefer Pathway~B for routine campaigns and manuscript figures.

---

## Shared thermodynamics

Both pathways use the same thermodynamic distance

$$
f=\alpha\frac{\Gamma^{\beta}}{1+\gamma\,r_s^{-\delta}\,\tau}\,,
$$

with production locks $\alpha=0.05$, $\beta=\gamma=\delta=1$. Pathway~A sets the polylog order $s=f$; Pathway~B uses the zeta argument $1+f$ in $W_\zeta=f\,\zeta(1+f)$. Weak coupling ($f\to 0$) returns ordinary RPA on both pathways.
