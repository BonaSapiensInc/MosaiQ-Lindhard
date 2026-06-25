# MosaiQ-Lindhard

Deterministic C++20 sandbox and companion tooling for finite-temperature Lindhard linear response and two-component RPA dielectric functions in homogeneous electron–ion plasmas.

Companion manuscript: *Linear response representation of two-fermion plasmas* (`manuscript/two-fermi.tex`).

## Repository layout

| Path | Purpose |
|------|---------|
| `simulator/` | MosaiQ-Lindhard C++20 engine (`mosaiq_simulator` CLI, CTest suite) |
| `scripts/` | Python figure generators (matplotlib / seaborn) |
| `manuscript/` | LaTeX source and cover letter |
| `docs/` | Architecture blueprint and legacy reference extracts |
| `output/` | CLI and plotting artifacts (gitignored; regenerate locally) |

## Quick start

### Build the simulator

```bash
cmake -S simulator -B simulator/build -DCMAKE_BUILD_TYPE=Release
cmake --build simulator/build
```

Run the CLI from the repository root:

```bash
./simulator/build/mosaiq_simulator
```

### Run tests

```bash
ctest --test-dir simulator/build --output-on-failure
```

### Regenerate figures

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

# Dynamic structure factor panels (requires output/output_structure_factor.dat)
python3 scripts/plot_Sqw.py

# Static S(q) gamma sweep (requires simulator run below)
./scripts/regenerate_sq_gamma.sh 1 10,50,100,150

python3 scripts/plot_plasmon_dispersion.py
python3 scripts/plot_t0_error.py
python3 scripts/plot_lindhard_t0_2d.py
```

Run the gamma sweep from the repository root (writes `output/output_Sq_gamma.dat`):

```bash
./simulator/build/mosaiq_simulator --gamma-sweep 1 10,50,100,150
python3 scripts/plot_Sq_gamma.py
```

### Compile the manuscript

From `manuscript/`, after all `output/*.pdf` figures exist:

```bash
cd manuscript && pdflatex two-fermi.tex && pdflatex two-fermi.tex
```

If `pdflatex` reports missing `S_*_contour.pdf` or `S_*_3d.pdf`, run `python3 scripts/plot_Sqw.py` first.

## Data availability

Source code, test suites, and build instructions are published under the **MosaiQ-Lindhard** name to avoid confusion with unrelated third-party packages (including the Lie-algebra library *Sophus*).

Repository URL (update before release):

`https://github.com/YOUR_GITHUB_USERNAME/MosaiQ-Lindhard`

## Citation

If you use this code, please cite the accompanying *Physical Review E* manuscript (in preparation).
