<p align="center">
  <!-- Absolute raw URL, not a relative path: this file is also the PyPI
       long_description, where relative links do not resolve. -->
  <img src="https://raw.githubusercontent.com/mbustama/oroscope/main/docs/source/_static/oroscope_logo.png"
       alt="Oroscope" width="200">
</p>

[![tests](https://github.com/mbustama/oroscope/actions/workflows/tests.yml/badge.svg)](https://github.com/mbustama/oroscope/actions/workflows/tests.yml)
[![Code Quality](https://github.com/mbustama/oroscope/actions/workflows/lint.yml/badge.svg)](https://github.com/mbustama/oroscope/actions/workflows/lint.yml)
[![codecov](https://codecov.io/gh/mbustama/oroscope/branch/main/graph/badge.svg)](https://codecov.io/gh/mbustama/oroscope)
[![Documentation](https://img.shields.io/badge/docs-GitHub%20Pages-blue.svg)](https://mbustama.github.io/oroscope/)
[![PyPI](https://img.shields.io/pypi/v/oroscope.svg)](https://pypi.org/project/oroscope/)
[![Downloads](https://pepy.tech/badge/oroscope)](https://pepy.tech/project/oroscope)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![Code style: ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

# Oroscope

**Terrain site search for particle-astrophysics observatories.** Greek *oros*, mountain,
and *skopein*, to look at.

> Coverage is not uploaded yet, so that badge stays grey until it is. The rest are live.

Oroscope searches digital elevation models for ground that can host an observatory. It
answers one structural question — *from this patch of ground, is there a target surface
at the right range, in the right direction, at the right relative orientation, with the
right matter behind it?* — which is what lets a single engine serve experiments that
look nothing alike. GRAND wants terrain a few degrees below the horizon, tens of
kilometres away; TAMBO wants a canyon wall two to five kilometres across. They differ in
their numbers, not their structure.

**Documentation:** [the physics](https://mbustama.github.io/oroscope/physics.html) ·
[assumptions and limitations](https://mbustama.github.io/oroscope/assumptions.html) ·
[quickstart](https://mbustama.github.io/oroscope/quickstart.html)

```bash
pip install oroscope
```

```python
from oroscope import site_searcher as ss

results = ss.run_from_config("config/grand_colca_config.json",
                             run_output_dir="output/colca")

print(results["explanation"])       # what was found, and why
```

or, the same search from a shell:

```bash
oroscope --config_path config/grand_colca_config.json
```

---

## Overview

Oroscope screens digital elevation models for ground that can host a
particle-astrophysics observatory, and reports where each experiment is viable, where
several are, and how much of the answer rests on assumptions.

It was written for **GRAND** (radio detection of air showers from Earth-skimming tau
neutrinos) and now serves **TAMBO** (particle detection across a deep canyon) through
the same scan engine, because the two ask the same structural question and differ in
their numbers. Adding an experiment means writing a configuration, not a code path.

Searching a DEM of hundreds of millions of pixels needs care with both memory and time,
so the tool uses **out-of-core memory mapping** throughout and a **Numba-compiled,
parallel** scan kernel. Terrain is screened by slope, aspect, altitude and exclusion
zones; the survivors are scanned over arrival directions; the results are scored against
per-experiment criteria and turned into sites with detector capacity.

---

## 1. Install and get a DEM

Python 3.9+.

```bash
pip install oroscope          # or, from a clone: pip install -e .
```

That pulls in the dependencies and installs six console scripts: `oroscope`,
`oroscope-combine`, `oroscope-crop`, `oroscope-sensitivity`, `oroscope-fetch-dem` and
`oroscope-fetch-roads`.

A search needs a **Digital Elevation Model** in `.tif` form, ~30 m resolution. Four
regions are bundled and fetch themselves, given a free
[OpenTopography](https://portal.opentopography.org/myopentopo) key:

```bash
oroscope-fetch-dem --region arequipa --open_topography_api_key YOUR_KEY
```

It writes the `.tif` into `input/dem/` and a matching configuration into `config/`,
both relative to the working directory unless `--output_dir` and `--config_dir` say
otherwise.

**[Installation and data →](https://mbustama.github.io/oroscope/installation.html)** for
conda environments, `oroscope.generate_env`, merging tiles for a region of your own, and
what each bundled DEM covers.

---

## 2. Quick-Start Guide

**Everything below can be done from Python, and that is the recommended way in.** The
command line is a thin wrapper over the same functions — argument parsing and file
placement, nothing else. There is no CLI-only behaviour, and anything the command line
can reach that the library cannot is a bug.

### Run a search

```python
from oroscope import site_searcher as ss

results = ss.find_grand_regions_interactive(
    dem_path="input/dem/colca.tif",
    run_output_dir="output/colca",
    search_mode="distributed", grid_type="hex",
    min_slope_deg=3.0, max_slope_deg=25.0,       # deployable ground
    min_dist_km=10.0, max_dist_km=40.0,          # where a tau may exit
    elev_min_deg=-3.0, elev_max_deg=3.0,         # the arrival window
    antenna_spacing_km=1.0, target_antennas=10000,
    downsample_factor=4, num_cores=8,
)
```

It **returns its results**, so nothing has to re-read the file it just wrote:

```python
print(results["results"]["total_sites"])       # how many sites
print(results["results"]["total_capacity"])    # how many detectors
print(results["explanation"])                  # the run, in plain language
```

### Start from a configuration rather than a signature

Configurations are data, and the template names every knob the tool understands:

```python
config = ss.default_config("arequipa")     # every key, with its default
ss.generate_config("arequipa.json", "arequipa")   # write it out
config = ss.load_config("arequipa.json")          # read one back

config["min_slope_deg"] = 5.0
results = ss.find_grand_regions_interactive(
    run_output_dir="output/arequipa",
    **{k: v for k, v in config.items()
       if not k.startswith("_")
       and k not in ("print_info", "output_directory_base_with_given_json")})
```

> **The three sources of defaults agree.** A parameter's default is the same whether you
> read it off the function signature, off `oroscope --help`, or out of
> `default_config()` — they disagreed on ten of them once, and a test now pins all three
> together. Starting from `default_config()` is still the clearer habit, because it puts
> every knob in front of you.

### Read what came back

```python
from oroscope import explain

chosen, shortlisted = explain.selected_sites(results)   # `sites` can exceed the selection
for site in chosen:
    print(site["site_id"], site["area_km2"], site["capacity_exact"],
          site["center_lat"], site["center_lon"])

binding = explain.binding_constraint(results["funnel"])
print(f"{binding['stage']} kept {100 * binding['kept_fraction']:.1f}%")
print(f"change: {binding['knob']}")

for entry in explain.site_strengths(chosen[0]["arrival_scan"]):
    print(entry["label"], entry["score"], entry.get("evidence"))
```

`explain.explain_results(results)` is a **pure function of the results dictionary** — no
DEM, nothing re-run — so a search from months ago can still be explained from its JSON:

```python
import json
with open("output/colca/oroscope_results_colca.json") as f:
    print(explain.explain_results(json.load(f)))
```

### Cut a DEM, combine experiments, test an assumption

```python
from oroscope import crop_dem, combine_experiments as combine, sensitivity

info = crop_dem.crop("input/dem/arequipa_SRTMGL1.tif", "input/dem/colca.tif",
                     north=-15.30, south=-15.85, west=-72.40, east=-71.55)

grand = combine.load_run("output/grand_colca_config")
tambo = combine.load_run("output/tambo_colca_config")
combine.check_alignment([grand, tambo])        # refuses to overlay the wrong ground

point = sensitivity.run_once(config, "output/sweep_point")
print(sensitivity.summarise(point))
```

### Guard the memory before a long run

```python
report = ss.preflight_memory("input/dem/arequipa_SRTMGL1.tif",
                             downsample_factor=4, candidate_stride=5)
print(report["estimate_gb"], report["available_gb"])
```

The estimate is what decides `downsample_factor` and `candidate_stride`: the full
Arequipa DEM is estimated at 7.2 GiB at 1 and 5.1 GiB at 4, against a measured **6.59 GiB
resident and 7.80 GiB of address space** at 4 — so read the estimate as a floor, and see
`--max_memory_gb` below before sizing a cap from it. Downsampling helps less than it looks —
it scales the labelling arrays as its inverse square but not the candidates, which are
taken on the native grid and dominate at this scale, so `candidate_stride` is the lever
on the larger term. Passing
`max_memory_gb` to a search caps its address space, so one that outgrows the machine
fails with `MemoryError` naming itself rather than inviting the OOM killer to pick a
victim.

### The same things from the command line

```bash
oroscope --generate_config arequipa.json --config_preset arequipa
oroscope --config_path arequipa.json
oroscope --config_path arequipa.json --min_slope_deg 5   # a typed flag beats the file
oroscope --config_path arequipa.json --resume --resume_dir output/arequipa
```

```bash
oroscope-crop input/dem/arequipa_SRTMGL1.tif input/dem/colca.tif \
    --north -15.30 --south -15.85 --west -72.40 --east -71.55
oroscope-combine output/grand_colca_config output/tambo_colca_config \
    --labels GRAND TAMBO --out output/combined
oroscope-sensitivity config/tambo_colca_config.json --sweep min_score 0.0 0.2 0.35 0.5
```

`origin_lat`/`origin_lon` are optional either way: the DEM's own tiepoint is used when
they are omitted, and a supplied value that disagrees with it is reported rather than
silently honoured.

### Output Products

By default, all generated output files are saved into a unified run folder located under `../output/`. If you use a JSON config file, the folder is named after the config file. If not, a timestamped folder (e.g., `../output/YYYYMMDD_HHMMSS/`) is automatically generated.

A complete run will produce the following files inside that directory:

* **`log.txt`**: A full transcript of the terminal execution, including settings used, memory usage, and runtime.
* **`explanation.txt`**: The run explained in plain language — what was found, which constraint set the size of the answer, what held each site back, and which numbers are assumptions. Written unless `--no_explain` is given.
* **`provenance.json`**: What produced the numbers — git commit and working-tree state, the DEM's sha256 and geometry, package versions, and the exact command.
* **`*.json`**: The results. Every resolved parameter, the **selection funnel** (survivors after each filter), the region accounting, per-stage timings, and a record for each site: area, capacity, facing direction, whether it was `selected`, and 34 aggregated scan observables plus each named score component. Note that `sites` lists everything clearing the thresholds while `total_sites` counts the selection — filter on `selected` before totalling.
* **`*.png`**: A high-resolution, annotated visualization map displaying the target terrain, overlaid RFI exclusion zones, and color-coded valid array sites. *(Format can be changed to PDF/SVG via parameters).*
* **`*.tif`**: A binary raster mask where `1` represents valid antenna deployment pixels and `0` is excluded terrain.
* **`*.tfw`**: An ESRI World File ensuring the `.tif` mask is properly georeferenced when loaded into GIS software (like QGIS or ArcGIS).
* **`*.kml`**: *(If `--generate_kml` is flagged)* A Google Earth compatible file containing the bounding polygons of all valid sites.

---

## 3. Development

The suite is standard-library `unittest` only, so it runs anywhere the tool does.
Terrain fixtures are synthetic with analytically known slope, aspect, target distance
and canyon geometry, so assertions are against arithmetic rather than a previous run.

```bash
cd tests && python -m unittest discover
```

Tests needing a real DEM skip when `input/dem/` is absent. After an intended change to
results, regenerate the golden files with `UPDATE_GOLDEN=1`.

**[Implementation notes →](https://mbustama.github.io/oroscope/implementation.html)**
for the benchmark harness, coverage of the Numba kernels, and the failure modes worth
knowing about. [docs/ROADMAP.md](https://github.com/mbustama/oroscope/blob/main/docs/ROADMAP.md) carries the development plan and every
measurement behind the current criteria.

---

## 4. Parameter Configuration Hierarchy

Parameters can come from four places. They are resolved in this order, **first match
wins**:

1. **An option you actually typed on the command line.** This beats everything, and says
   so when it overrides a config file. It used to *lose* to the config file, silently —
   and since `--generate_config` writes every key, a generated config made every flag on
   the command line a no-op with no warning.
2. **The config file** given with `--config_path`.
3. **`../config/fallbacks.json`**, if present. Useful for lab-wide defaults. Every value
   taken from here is announced, because a fallback is the least visible input the tool
   has.
4. **The built-in default**, as listed below.

Nothing here is required except a DEM: `origin_lat`/`origin_lon` are read from the DEM's
own GeoTIFF tiepoint when omitted, and a supplied origin that disagrees with the file by
more than ~100 m is reported rather than silently honoured.

### The options

There are 87 of them, and the complete reference — every flag with its type, default and
what it does — lives on the **[CLI page](https://mbustama.github.io/oroscope/cli.html)**.
It is generated against the parser and a test fails if the two drift, which a copy here
could not promise: this README carried its own table until it had quietly fallen ten
flags behind and was documenting one that no longer existed.

The four worth knowing before the first run:

| Option | What it decides |
| --- | --- |
| `--config_path` | The configuration to run. Everything else has a default; an explicitly typed option beats the file. |
| `--candidate_stride` | The memory and time lever. Unbiased in acceptance, but it costs area unless the closing element outruns the gap it leaves. |
| `--downsample_factor` | The resolution area and sites are measured at. Costs a thin feature more than a blocky one. |
| `--max_memory_gb` | An address-space cap, and **not the same quantity the estimate reports.** The estimate is anonymous memory; this bounds virtual address space, which on the full Arequipa DEM is 7.80 GiB against a 6.59 GiB resident peak. Sizing the cap from the estimate is how a run dies 25 minutes in. The default is 80% of what is free, which is below what that run needs. |

## 5. Internal Workings: The 6-Step Pipeline

The script processes terrain logically through six distinct architectural phases.

### Step 1: Disk Setup & Memory Management

To handle massive DEM files (which can easily exceed 20GB of RAM if loaded natively), the script instantly converts the input `.tif` into a Numpy `.npy` file. It then uses `np.lib.format.open_memmap` to establish "Ping-Pong" buffers (`buffer_A.npy`, `buffer_B.npy`) on the hard drive. All subsequent operations read and write to the disk in chunks, allowing the script to run seamlessly on standard laptops.

The map resolution is read from the DEM's `ModelPixelScaleTag` at this stage and reported in the run banner. Every downstream conversion — slope gradients, ray-tracing step lengths, RFI radii, morphology kernels, grid packing, and the georeferencing of the `.tif`/`.tfw`/`.kml` products — derives from it, so DEMs at resolutions other than 1 arc-second are handled correctly. Pass `--cell_size_deg` to override a DEM whose metadata is missing or wrong.

A geographic (EPSG:4326) DEM has pixels that are square in **degrees**, not in metres: a degree of longitude shrinks as `cos(latitude)`, so a 1 arc-second pixel spans about 30.7 m north-south but only ~29.5 m east-west in southern Peru. The pipeline therefore carries two metric pixel sizes and applies each on its own axis; angular quantities (the world file, KML coordinates, map axes) use the single degree value. The longitude scale is evaluated once at the DEM's centre latitude, which spreads the residual error of ignoring its north-south variation evenly across the map — about ±0.7% over a 3° tall DEM.

### Step 2: Topographic Screening

The code steps through the DEM in defined RAM chunks (configured by `--tile_size`).

1. Uses `np.gradient` to establish raw `dy` and `dx` vectors.
2. Derives physical `slope` and `aspect` angles using trig arrays.
3. Filters the terrain by the bounds (`min_slope`, `altitude`, `aspect`, etc.).
4. Evaluates geographic spatial logic. RFI exclusion zones are tested by real ground distance in metres, using the separate north-south and east-west pixel sizes, so a zone stays a true circle on the ground rather than becoming an ellipse.
5. Surviving pixels are thinned by `--candidate_stride` (default 5x) and passed forward as raw candidate coordinates.

### Step 3: The arrival scan

The expensive step, and the heart of the tool. For every surviving pixel it traces rays
*backwards* along a fan of arrival directions and asks what each one meets.

* **One walk per (candidate, azimuth).** Writing the terrain's elevation angle at ground
  distance `d` as `atan((z(d) - d²/2R - z₀)/d)`, a ray at angle θ first meets terrain at
  the smallest `d` where that exceeds θ. Because the running maximum only increases,
  each new maximum claims a contiguous band of elevation bins — so a single pass fills
  every bin at once. Elevation binning is therefore nearly free, and the **azimuth count
  is what sets the cost**.
* **Column depth from the same samples.** The ray is underground wherever the terrain
  angle exceeds θ, so binning the terrain angle and taking a suffix sum gives the
  underground path length for every bin, accumulating all the rock a ray crosses rather
  than only the first chord.
* **Compiled and parallel.** The kernel is Numba-compiled and spread across cores with
  `prange`, with candidates dealt in blocks so threads get comparable work without
  losing memory locality.
* **Two Earth radii.** Particles are not refracted, so the geometry uses the true
  6371 km; the radio path uses the 4/3 convention, and only for the Fresnel term.

The scan reports per-candidate observables — accepted solid angle, distance to the exit
point, column depth, horizon, atmospheric depth, Earth chord, far-wall slope — which are
then scored. See [the physics](https://mbustama.github.io/oroscope/physics.html) for the
derivation of each criterion.

### Step 4: Spatial Pruning

A single pixel that sees a mountain is useless if a truck cannot deploy an antenna there. The script uses SciPy's morphological kernels (`binary_closing`, `binary_opening`) on the massive memory maps.

* **Closing:** Fills in small, unviable gaps (potholes) in otherwise good slopes.
* **Opening:** Erases narrow, thin ridge lines (tendrils) that do not meet the `min_width_km` requirement.

### Step 5: Capacity Analysis

The script isolates disconnected sub-arrays using `scipy.ndimage.label`.
Instead of estimating area, it invokes `count_grid_capacity` which physically simulates dropping bounding boxes in either a staggered `hex` pattern or a strict `square` grid. Sites that cannot fit the `min_sub_array_size` are discarded. Capacity is counted **within each region**, not over its bounding box — which also contains other sites, and inflated the count by 38% on a canyon network.

Surviving sites are ranked by capacity and selected in that order. With
`--stop_at_target` selection stops once `target_antennas` is met, so the run reports the
best sites for the array actually wanted rather than every patch of qualifying ground.
Sites that qualified but were not selected stay in the results file, flagged
`selected: false`.

### Step 6: Output Generation

Everything is exported into a unified, dynamically generated run directory.

* **GeoTiff:** The final binary mask is saved alongside a `.tfw` (World File), allowing direct drag-and-drop into QGIS/ArcGIS.
* **KML:** Bounding polygons are extracted using Matplotlib's contour tool and formatted as yellow overlays for Google Earth.
* **JSON:** Every resolved parameter, the selection funnel, region accounting, per-stage timings, and a full record for each site.
* **Provenance:** Git commit, DEM checksum, package versions and the command, in their own file so they stay readable beside the science outputs.
* **Explanation:** The run in plain language, printed and saved. See §2's output list.