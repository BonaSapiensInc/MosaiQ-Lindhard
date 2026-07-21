# MosaiQ-Lindhard Documentation Map

Documentation is organized by **engine phase**. Phase 1 is the live two-fermion plasma / dual-pathway stack. Phase 2 is the electron–phonon Berezinian (superdeterminant) extension — documented before C++ scaffolding.

---

## Phase 1 — Two-fermion plasma (live)

| Document | Role |
|----------|------|
| [`theory.md`](theory.md) | Dual-pathway strong-coupling regularization (A = PolyLog-RPA, B = LFC–Zeta) |
| [`usage.md`](usage.md) | CLI pathways and invocation examples |
| [`simulator_architecture.md`](simulator_architecture.md) | C++20 engine blueprint (`core` / `physics` / `engine`) |
| [`zeta_rpa_integration.md`](zeta_rpa_integration.md) | Zeta-RPA insertion doctrine and Phase Z1–Z4 status |
| [`release_notes_v2.md`](release_notes_v2.md) | `v2.0.0-dual-pathway` release notes |

**Reading order (Phase 1):** `theory.md` → `usage.md` → `simulator_architecture.md` → `zeta_rpa_integration.md`.

---

## Phase 2 — Electron–phonon Berezinian (blueprint)

| Document | Role |
|----------|------|
| [`BEREZINIAN_ARCHITECTURE.md`](BEREZINIAN_ARCHITECTURE.md) | Graded supergraph / Schur / $DX{=}C$ / $i\eta$ / `std::span` numerical doctrine |

**Planned headers (not yet in tree):** `core/supermatrix.hpp`, `math/linalg_utils.hpp`, `lindhard/berezinian.hpp` — see §4 of the Berezinian architecture note. Mapping onto the live `simulator/src/{core,physics,engine}/` layout is deferred until scaffolding is authorized.

**Doctrine continuity:** Phase 2 inherits Phase 1 causality (real-axis / Kramers–Kronig / sinc) and extends the matrix response to a fermion–boson graded supermatrix. Zeta-RPA dresses the electron ladder; the Berezinian controls electron–phonon topology.

---

## Supporting material

| Path | Role |
|------|------|
| [`legacy_extraction/`](legacy_extraction/) | Historical C++ extracts (reference only; not the live build) |
| [`references/`](references/) | External figures and PDFs used during theory development |

---

## Manuscript

Companion paper: [`manuscript/two-fermi.tex`](../manuscript/two-fermi.tex). Architecture docs must not contradict the frozen dual-pathway numbering in Sec.~V / Appendix.
