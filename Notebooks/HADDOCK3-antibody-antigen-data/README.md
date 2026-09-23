# HADDOCK3 antibody-antigen tutorial data

Data for `../HADDOCK3-antibody-antigen-lite.ipynb`, the in-browser (JupyterLite)
version of the [HADDOCK3 antibody-antigen
tutorial](https://www.bonvinlab.org/education/HADDOCK3/HADDOCK3-antibody-antigen/).

This is a **trimmed copy** (1.6 MB) of the full tutorial bundle
(`HADDOCK3-antibody-antigen-notebook.zip`, 54 MB, linked from the Colab
notebook), containing only the files the notebook reads. It is shipped in the
repository because the browser cannot download the original archive: the
surfdrive host sends no `Access-Control-Allow-Origin` header, so the fetch is
blocked by the browser's cross-origin policy.

```
pdbs/          4G6K (antibody) and 4I1B (antigen) prepared for docking,
               the 4G6M reference complex and its energy-minimised S150W mutant
restraints/    active/passive residue lists and the AIR + body restraint files
workflows/     the three haddock3 configuration files used in the tutorial
scripts/       the analysis shell scripts (re-implemented in Python in the notebook)
runs/run1/                     pre-calculated docking run (caprieval statistics,
                               top cluster models, contact map, html report)
runs/run-energetics-alascan/   pre-calculated alanine scanning (BONUS 1)
runs/run-scoring/              pre-calculated scoring run (BONUS 2)
```

The notebook never writes here: anything it produces goes to
`../HADDOCK3-antibody-antigen-output/`.
