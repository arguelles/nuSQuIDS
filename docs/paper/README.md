# nuSQuIDS-GPU paper — draft

`paper.pdf` is a WIP draft describing the CUDA GPU backend added to nuSQuIDS
on the `cuda-backend` branch. It is included here so PR reviewers have
access to the full validation and performance story alongside the code.

## Contents

- `paper.pdf` — the compiled draft (RevTeX 4-2, PRD reprint style)
- `paper.tex` + `sections/*.tex` — LaTeX sources
- `figures/*.pdf` — plots: `speedup.pdf`, `residuals.pdf`, `perbin_tau*.pdf`,
  `perbin_icecube*.pdf`
- `data/*.csv` — benchmark measurements (CPU/GPU wall-time by grid size on
  MIG 1g.10gb and full A100-SXM4-80GB), interaction-mode residuals
- `scripts/*.py` — plotting scripts that convert `data/*.csv` into `figures/*.pdf`

## Highlights from the draft

- Validation: five representative scenarios agree with the CPU reference at
  the $\sim 10^{-3}$ relative level; per-bin output is byte-identical between
  two A100 partitions of the same cluster.
- Performance: up to $48\times$ speedup at $(n_{c_\theta}, n_E) = (200, 200)$
  with full non-coherent interactions (NC + CC + tau regen + Glashow) on a
  full A100-SXM4-80GB — $2\times$ better than the previously published GPU
  port (Kallenborn et al. 2019).
- Known limits documented as future work: multi-path-per-block scheduling
  and an embedded fifth-order integrator to close the gap on
  low-$n_{c_\theta}$ IceCube-shape workloads.

## Rebuilding

```sh
cd docs/paper
latexmk -pdf paper.tex
```

`references.bib` and `sections/*.tex` are pulled in automatically.

## Regenerating figures

The Python scripts require `matplotlib`, `numpy`, and `pandas`:

```sh
cd docs/paper
python scripts/plot_speedup.py       # figures/speedup.pdf
python scripts/plot_residuals.py     # figures/residuals.pdf
python scripts/plot_perbin.py        # figures/perbin_tau*.pdf, figures/perbin_icecube*.pdf
```

Each script reads from `data/` and writes to `figures/`.

## Status

Draft. The A100 sweep is complete; H100 / H200 and multi-GPU strong-scaling
rows are still to be measured. Real-world atmospheric-analysis example is
still to be added. The `todo` macros in `paper.tex` mark open items.
