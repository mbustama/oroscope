# Stored results — the full Arequipa DEM

This directory holds the *small* artefacts of the three full-DEM searches, so that
[notebook 7](../../notebooks/07_running_a_search.ipynb) can **read** them instead of
running them.

That split is the point. Each of these searches takes about half an hour, and CI
executes every notebook on every push; a tutorial that costs ninety minutes of compute
per commit is a bill rather than a tutorial. So the expensive half runs once, locally,
on a machine that has the DEM — and what lands here is a few hundred kilobytes of JSON
that opens instantly and needs no DEM at all. That is possible because
`explain.explain_results()` is a pure function of the results dictionary.

## Producing it

```bash
python tools/run_arequipa_full.py --dry-run   # what it will do, and what it will cost
python tools/run_arequipa_full.py             # GRAND, TAMBO, then the combination
```

**Regenerate when a configuration changes, and not otherwise.** `manifest.json` records
when the store was built, from which DEM and which configs; each `*_provenance.json`
records the git commit, the DEM's sha256 and the package versions. A stale store is
therefore detectable rather than merely suspected — which matters, because the whole
premise of storing rather than recomputing is that nobody looks again.

## What lands here

| file | what |
| --- | --- |
| `manifest.json` | When, from what, by what. Read this first. |
| `grand_results.json` | GRAND's full results: parameters, funnel, regions, timings, per-site records. |
| `grand_provenance.json` | Commit, DEM checksum, package versions, command. |
| `grand_explanation.txt` | The run in plain language. |
| `tambo_*` | The same three, for TAMBO. |
| `combined_report.json` | Joint, union and co-location over the two masks. |

**Not** the rasters. A GeoTIFF mask of a 129 Mpx DEM is far too large for a repository,
and the notebook does not need one. The full outputs — GeoTIFF, world file, KML, PNG,
log — stay in `output/`, which is gitignored.

## The configurations

`config/grand_arequipa_full.json` and `config/tambo_arequipa_full.json`. They are the
Colca crop configs with three changes: the full DEM, the origin read from the file's own
tiepoint rather than the crop's corner, and `downsample_factor: 4` instead of 1 — the
estimator puts this DEM at 4.5 GiB of anonymous memory at 1 against ~6 GiB typically
free, and 2.3 GiB at 4.

Every criterion is otherwise unchanged from the crop, deliberately: the point of this
run is **scale**, not a different question. A crop is chosen because it is interesting,
and a search over ground chosen for being interesting is not a survey.

That `downsample_factor` has a price worth stating where the numbers are read: area is
measured on the downsampled mask while capacity is measured at full resolution, so a
feature a few pixels wide keeps its detectors and loses area. Read these areas as lower
bounds, and more so for TAMBO's canyon strips than for GRAND's blobs.
