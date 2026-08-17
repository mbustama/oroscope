# Stored results — the full Arequipa DEM

This directory holds the *small* artefacts of the three full-DEM searches, so that
[notebook 9](../../notebooks/09_arequipa_dem.ipynb) can **read** them instead of
running them.

That split is the point. Each of these searches takes about half an hour, and CI
executes every notebook on every push; a tutorial that costs ninety minutes of compute
per commit is a bill rather than a tutorial. So the expensive half runs once, locally,
on a machine that has the DEM — and what lands here is a few hundred kilobytes of JSON
that opens instantly and needs no DEM at all. That is possible because
`explain.explain_results()` is a pure function of the results dictionary.

## Producing it

```bash
python tools/run_full_dem.py --region arequipa --dry-run              # the cost, then stop
python tools/run_full_dem.py --region arequipa --max-memory-gb 8.0    # both, then combine
python tools/run_full_dem.py --region arequipa --max-memory-gb 8.0 --only grand
```

The cap is spelled out rather than left to the default for the reason given under
[The configurations](#the-configurations): 80% of available lands *below* what this
search needs.

`--dry-run` **starts nothing** — no search, no file, no change to this directory. It
prints the five things worth knowing before committing an hour of a machine:

```text
region:    arequipa
DEM:       input/dem/arequipa_SRTMGL1.tif              which file, and whether it is there
estimate:  5.56 GiB at downsample_factor 4, candidate_stride 5
           5.08 search + 0.48 map, then 1.32 to combine   where the estimate goes
available: 7.0 GiB                                     against what the system reports
would run: grand, tambo, then combine                  honouring --only
expected:  ~25 min for grand, ~1 min for tambo         so you do not start it blind
store:     results/arequipa_full                       where the artefacts land
```

The estimate is what decides `downsample_factor` and `candidate_stride`: the same DEM
needs 7.2 GiB at 1 and 5.1 GiB at 4. Downsampling scales the labelling arrays as its
inverse square but not the candidates, which are taken on the native grid and dominate
at this scale, so striding is the stronger lever. It deliberately
excludes the memory-mapped DEM, which is file-backed and evictable — counting it would
make every large search look impossible when the streaming design exists precisely so
that it is not. No memory cap is applied during a dry run, since nothing is allocated.

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
| `combination_explanation.txt` | The combination in plain language. Copied beside the report rather than left behind: when it was not, a re-combine refreshed the numbers while the prose beside them kept the old ones, and the region notebooks print this file verbatim. |

**Not** the rasters. A GeoTIFF mask of a 129 Mpx DEM is far too large for a repository,
and the notebook does not need one. The full outputs — GeoTIFF, world file, KML, PNG,
log — stay in `output/`, which is gitignored.

## The configurations

`config/grand_arequipa_full.json` and `config/tambo_arequipa_full.json`. They are the
Colca crop configs with three changes: the full DEM, the origin read from the file's own
tiepoint rather than the crop's corner, and `downsample_factor: 4` instead of 1 — the
estimator puts this DEM at 7.2 GiB of anonymous memory at 1 against ~6-7 GiB typically
free, and 5.1 GiB at 4. Even at 4 it needs a machine whose desktop is not holding half
of RAM: the run measured **6.59 GiB peak RSS against 7.80 GiB of address space**, and
wants `--max-memory-gb` set explicitly rather than the default 80%-of-available cap.
Size it from the virtual figure — that is what `RLIMIT_AS` bounds. A 6.5 GiB cap, above
the resident peak and below the virtual one, still died in the final analysis.

Every criterion is otherwise unchanged from the crop, deliberately: the point of this
run is **scale**, not a different question. A crop is chosen because it is interesting,
and a search over ground chosen for being interesting is not a survey.

That `downsample_factor` has a price worth stating where the numbers are read: area is
measured on the downsampled mask while capacity is measured at full resolution, so a
feature a few pixels wide keeps its detectors and loses area. Read these areas as lower
bounds, and more so for TAMBO's canyon strips than for GRAND's blobs.
