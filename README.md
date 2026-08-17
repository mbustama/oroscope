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

> Several of these badges are not live yet: coverage is not uploaded, the package is not
> on PyPI, and the documentation is not deployed to Pages. They are in place so that
> turning each of those on needs no edit here.

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
from oroscope import site_searcher as ss, explain

results = ss.find_grand_regions_interactive(
    dem_path="input/dem/colca.tif", run_output_dir="output/colca",
    **ss.load_config("config/grand_colca_config.json"))

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

## 1. Requirements

The package requires **Python 3.9+**. Due to the heavy reliance on C-compiled math and geospatial array processing, using a virtual environment (like Conda) is highly recommended.

### Core Dependencies:

* `numpy` (Core matrix math)
* `scipy` (Morphological image processing and connected-component labeling)
* `numba` (JIT compiler for parallelized physics kernels)
* `tifffile` (Reading/Writing large geospatial TIFFs)
* `matplotlib` (Generating visualization maps and extracting KML contours)
* `tqdm` (Progress bars)
* `psutil` (Optional, used for printing system RAM diagnostics)

**Installation via pip:**

```bash
pip install oroscope          # or, from a clone: pip install -e .
```

That installs the dependencies and five console scripts: `oroscope`,
`oroscope-combine`, `oroscope-crop`, `oroscope-sensitivity` and `oroscope-fetch-dem`.

### Automated Conda Environment Generation:

Instead of installing packages manually, you can use the included `generate_env.py` script. This tool parses the main script using Python's Abstract Syntax Tree (AST) to securely identify all required third-party dependencies, checks what is missing from your active environment, and generates a clean, `conda-forge` prioritized `environment.yml` file.

**Usage:**

```bash
# Generate environment.yml from what the code actually imports.
# Run from src/: it reads site_searcher.py in the working directory.
cd src && python generate_env.py

# Create the new Conda environment from the generated file
conda env create -f environment.yml

# Activate the environment
conda activate oroscope

```

It parses the sources rather than importing them, so it works in an environment where
the dependencies are precisely what is not installed yet. `pyproject.toml` remains the
authoritative list.

---

## 2. Setup & Data Acquisition

The script requires a **Digital Elevation Model (DEM)** in `.tif` format, optimized for ~30-meter resolution models. We highly recommend using **[OpenTopography](https://opentopography.org/)** to acquire this data (e.g., SRTM or ALOS AW3D30).

### Automated Setup (Recommended)

We provide a `fetch_dem.py` script that verifies your environment dependencies, automatically downloads the required DEM files for the primary target regions (Lima and Arequipa), and generates ready-to-use configuration files.

**Step 1: Obtain an OpenTopography API Key**

1. Create a free account at [OpenTopography](https://portal.opentopography.org/myopentopo).
2. Log in and navigate to the **"myOpenTopo"** dashboard.
3. Click on **"Request an API Key"** to generate your unique authorization token.

**Step 2: Run the Setup Script**
Pass your API key to the setup script to begin the automated download and configuration process:

```bash
oroscope-fetch-dem --open_topography_api_key YOUR_API_KEY_HERE

```

This downloads the `.tif` files into `../input/dem/` and generates matching JSON config
files in `../config/`. Run it from `src/`, since both paths are relative to it.

### Manual Setup

If you are targeting a region other than Lima or Arequipa:

* Download the required regional tiles manually via the OpenTopography web portal.
* **Preparation:** If your target region spans multiple tiles, merge them into a single `.tif` using a GIS tool like QGIS or GDAL (`gdal_merge.py`) before running the script.

---

## 3. Quick-Start Guide

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
Arequipa DEM needs 7.2 GiB at 1 and 5.1 GiB at 4. Downsampling helps less than it looks —
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

## 3b. Development: tests and benchmarks

The test suite uses only the standard library's `unittest`, so it runs anywhere the
tool itself runs — no extra dependencies. Terrain fixtures are synthetic with
analytically known slope, aspect, target distance and canyon geometry, so assertions
are against arithmetic rather than against a previous run.

```bash
cd tests && python -m unittest discover
```

Tests that need a real DEM skip automatically when `input/dem/` is absent (it is
gitignored). After an intended change to results, regenerate the golden files:

```bash
cd tests && UPDATE_GOLDEN=1 python -m unittest test_regression
```

The benchmark harness records per-stage wall time and peak memory on fixed inputs and
compares against `bench/baseline.json`, failing if any stage slows by more than 30%:

```bash
python bench/benchmark.py
```

Every run also writes a **selection funnel** to the log and results JSON, showing how
many pixels survived each filter, plus a `provenance.json` capturing the git commit,
DEM checksum, resolved parameters and package versions. When a search returns no
sites, the funnel is the first place to look.

And every run **explains itself**, in plain language, at the end and in
`explanation.txt`: what was found, which funnel stage set the size of the answer and
which parameter to change, which named score component held each site back, and which
of the numbers are assumptions rather than measurements. On by default; `--no_explain`
suppresses it. Any results file can be re-explained later without re-running anything:

```python
import json, explain
print(explain.explain_results(json.load(open("output/.../oroscope_results_colca.json"))))
```

See [docs/ROADMAP.md](docs/ROADMAP.md) for the development plan and a review of the
current selection criteria.

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

### Complete List of Options

Every option the tool accepts. Anything shaping a result is here — nothing that matters
is hard-coded.

#### Required Data Inputs

| Argument | Type | Default | Description |
| --- | --- | --- | --- |
| `--dem_path` | String | — | Path to the input elevation `.tif`. The only genuinely required input. |
| `--origin_lat` | Float | from DEM | Latitude of the DEM's north-west corner. Read from the GeoTIFF tiepoint when omitted, which is the recommended use. |
| `--origin_lon` | Float | from DEM | Longitude of the same corner. |

#### Physical Layout & Geometry

| Argument | Type | Default | Description |
| --- | --- | --- | --- |
| `--target_antennas` | Int | `10000` | Capacity wanted from the array. |
| `--antenna_spacing_km` | Float | `1.0` | Distance between detectors. |
| `--grid_type` | String | `hex` | Deployment lattice, `hex` or `square`. |
| `--min_width_km` | Float | `2.0` | Narrowest feature to keep. `0` disables pruning, which is what a strip along a canyon wall needs. |
| `--min_sub_array_size` | Int | `500` | Capacity a disconnected sub-array must reach to count. |
| `--search_mode` | String | `distributed` | `single` finds one monolithic site; `distributed` allows sub-arrays. |
| `--stop_at_target` | Flag | off | Stop selecting sites once `target_antennas` is met. Sites are ranked by capacity, so this reports the best sites for the array actually wanted rather than every patch of qualifying ground. |
| `--gap_close_km` | Float | `antenna_spacing_km` | Morphological closing element. Closing more than doubled the reported area on real terrain (2.29× at Colca), so this is worth setting deliberately; `0` disables it. |

#### Topographic Screening

| Argument | Type | Default | Description |
| --- | --- | --- | --- |
| `--min_slope_deg` | Float | `3.0` | Minimum terrain steepness. |
| `--max_slope_deg` | Float | `25.0` | Maximum terrain steepness. |
| `--slope_baseline_m` | Float | DEM resolution | Ground distance over which slope is measured. It matters: the same terrain gives a median 17.8° at ~61 m and 10.8° at 1 km. |
| `--min_altitude` | Float | `None` | Minimum site altitude, in metres. |
| `--max_altitude` | Float | `None` | Maximum site altitude, in metres. |
| `--min_aspect_deg` | Float | `None` | Restrict sites to a facing direction (0–360). |
| `--max_aspect_deg` | Float | `None` | Upper bound of that range. |
| `--candidate_stride` | Int | `5` | Keep every Nth screened pixel before scanning. Measured unbiased against a stride-1 control. Use `1` to scan every candidate. |

#### The Arrival Scan

| Argument | Type | Default | Description |
| --- | --- | --- | --- |
| `--min_dist_km` | Float | `10.0` | Nearest accepted exit point. |
| `--max_dist_km` | Float | `80.0` | Furthest accepted exit point. |
| `--energy_min_pev` | Float | `None` | Lower tau energy, in PeV. With `--energy_max_pev`, derives the decay-baseline distance window instead of setting it by hand. |
| `--energy_max_pev` | Float | `None` | Upper tau energy, in PeV. |
| `--elev_min_deg` | Float | `-3.0` | Lower edge of the accepted arrival-elevation window. |
| `--elev_max_deg` | Float | `3.0` | Upper edge of that window. |
| `--n_elev_bins` | Int | `12` | Bins across the window. Nearly free: one profile walk serves every bin. |
| `--n_azimuths` | Int | `9` | Azimuths per candidate. **This is what sets the cost.** |
| `--azimuth_half_width_deg` | Float | `60.0` | Half-width of the azimuth fan about the pixel's aspect. `-1` sweeps the full 360°. |
| `--max_range_km` | Float | `max_dist_km` | How far each profile is walked. Worth setting larger for a short-range search, or the reported column depth is a property of where the walk stopped rather than of the target. |
| `--min_column_depth_gcm2` | Float | `0.0` | Column depth a direction must have to count. |
| `--min_target_slope_deg` | Float | `None` | Require the struck terrain to be at least this steep, along the arrival azimuth. This is what separates a canyon *wall* from a hillside. |
| `--max_target_slope_deg` | Float | `None` | Upper bound on the struck terrain's slope. |
| `--muon_shielding_km` | Float | `None` | Rock overburden required to reject atmospheric muons (TAMBO quotes >4). A floor on column depth, not a band. |
| `--require_sky` | Flag | off | Invert the test: accept directions reaching clear sky, for cosmic-ray-style channels. |
| `--nearest_sampling` | Flag | off (bilinear on) | Sample profiles at pixel centres instead of interpolating. Faster, but treats terrain as blocky, which over-estimates blocking. Measured: bilinear gives +13.4% acceptance at 1.44× cost. |
| `--refraction_k` | Float | `4/3` | Refraction k-factor, **radio path only**. Particle trajectories always use the true Earth radius. |

#### Radio Propagation

| Argument | Type | Default | Description |
| --- | --- | --- | --- |
| `--fresnel_frequency_mhz` | Float | `None` | Band for the Fresnel clearance measurement, e.g. 50. Omitting it skips that pass entirely. |
| `--antenna_height_m` | Float | `2.0` | Antenna height above ground, for that measurement. |
| `--fresnel_near_field_m` | Float | `500.0` | Skip this much of the path when measuring clearance. Below ~500 m the measure is dominated by the ground beside the antenna. |
| `--include_near_field` | Flag | off (excluded) | Measure from the antenna outward instead. For study only, for the reason above. |
| `--geomag_declination_deg` | Float | `-6.9` | Declination, degrees east of north. **Does not follow the site** — supply the IGRF value per region. |
| `--geomag_inclination_deg` | Float | dipole | Inclination, positive downward. Defaults to a centred-dipole estimate at the DEM's own centre, so it does follow the site. |
| `--no_geomagnetic` | Flag | off (weighting on) | Ignore the geomagnetic angle and weight all directions equally. Radio emission goes as \|v × B\|; particles do not care. |

#### Shower and Scoring

| Argument | Type | Default | Description |
| --- | --- | --- | --- |
| `--min_score` | Float | `0.0` | Discard candidates scoring below this. **A product score has no safe threshold** — measured, 0.0/0.35/0.5 gave 45928/2056/0 detector positions. Prefer the percentile below. |
| `--score_percentile` | Float | `None` | Keep this percentage of viable candidates by rank instead. Scale-free, and preferred for exactly that reason. |
| `--score_composition` | String | `product` | How components combine: `product`, `mean` or `min`. |
| `--score_weights` | String | `None` | Per-component weights for `weighted` composition, as `shower=2,solid_angle=1`. |
| `--grammage_mode` | String | `radio` | How atmospheric depth is scored: `radio` is a maturity threshold (emission comes from shower maximum and then propagates through transparent air); `particle` is a band (particle content dies after maximum). |
| `--grammage_maturity_gcm2` | Float | `700` | Depth at which the `radio` ramp reaches 1. |
| `--grammage_band_gcm2` | Float ×2 | `(700, 2800)` | Depth band scoring 1 in `particle` mode. A short crossing gives far less — Colca supplies ~170 g/cm², so this must be lowered there or nothing scores. |
| `--grammage_band_fraction` | Float | `0.1` | Fraction of peak particle content still counting as a usable shower, when the band is derived from an energy range. |
| `--shower_elongation_rate_gcm2` | Float | `55` | How much deeper shower maximum sits per decade of energy. 85 for a purely electromagnetic cascade. |
| `--shower_lambda_gcm2` | Float | `70` | Gaisser–Hillas interaction length, setting how fast the profile rises and falls. |
| `--shower_development_m` | Float | `3000.0` | Path the shower needs after the tau decays. |
| `--depth_band_gcm2` | Float ×2 | `None` | Column-depth band scoring 1. A band, not a floor: the tau must be produced *and* escape. |
| `--distance_band_m` | Float ×2 | decay window | Exit-distance band scoring 1. |
| `--solid_angle_half_fraction` | Float | `0.076` | Fraction of the sky the azimuth fan and arrival window could accept that scores 0.5. Dimensionless, so it survives a change of fan width or arrival window untouched. |
| `--solid_angle_half_sr` | Float | `None` | The same thing as an absolute solid angle, in steradians, for a caller who wants one. Not portable: it silently encodes the fan width and the elevation window, which is why the fraction is the default. |
| `--clearance_full_at` | Float | `1.0` | Fresnel clearance ratio, in first-Fresnel radii, that scores 1. |
| `--nu_interaction_length_gcm2` | Float | `None` | Neutrino interaction length, enabling the Earth-chord attenuation term (order 1e8 near an EeV). |
| `--decay_energy_pev` | Float | `None` | Score the decay probability at a single tau energy. Left off by default: across one experiment's reach this *chooses* the answer rather than approximating it. |
| `--decay_energy_min_pev` | Float | `None` | Lower end of the energy range for the decay term. With the maximum, folds the probability over a power-law spectrum — the defensible form. |
| `--decay_energy_max_pev` | Float | `None` | Upper end of that range. |
| `--decay_spectral_index` | Float ×1–2 | `2.0` | Spectral index γ for the folded decay term. One value pins the spectrum; two marginalise uniformly over the range, which says "not known" rather than pretending to a value. |

#### Logistics & Exclusions

| Argument | Type | Default | Description |
| --- | --- | --- | --- |
| `--rfi_zones` | String | `none` | Radio exclusion zones: a preset (`lima`, `arequipa`) or a JSON list of polygons/circles. |
| `--road_map_path` | String | `None` | Aligned `.tif` of distance-to-roads. |
| `--max_road_dist_km` | Float | `20.0` | Maximum allowed distance from a road, when a road map is given. |

#### Compute & System Management

| Argument | Type | Default | Description |
| --- | --- | --- | --- |
| `--num_cores` | Int | `-1` | CPU threads for the scan. `-1` uses all available. Measured scaling: 1.85× at 2 threads, 3.70× at 8. |
| `--tile_size` | Int | `2048` | Square chunk loaded into RAM at a time. Reduce if memory is tight. |
| `--downsample_factor` | Int | `4` | Coarsening applied before labelling and area measurement. Memory scales as its inverse square — the knob that matters most for a large DEM. |
| `--cell_size_deg` | Float | from DEM | Map resolution in degrees per pixel. Read from the GeoTIFF's own tags; set only to override a DEM with wrong metadata. |
| `--max_memory_gb` | Float | 80% of free | Ceiling on the process's address space, in GiB, so a search that outgrows the machine fails with `MemoryError` rather than inviting the OOM killer to pick a victim. `0` disables the cap. |
| `--resume` | Flag | off | Resume from the ray-tracing checkpoint if buffers exist. |
| `--resume_dir` | String | run dir | Path to a failed run's directory to resume from. |

#### Output & Metadata

| Argument | Type | Default | Description |
| --- | --- | --- | --- |
| `--region_name` | String | `None` | Cosmetic name printed on the visualization map. |
| `--generate_kml` | Flag | off | Also write a Google Earth `.kml` of the findings. |
| `--output_image_format` | String | `png` | Format of the map: `png`, `pdf`, `svg`. |
| `--output_directory_base_with_given_json` | String | `../output/` | Base directory for run folders. |
| `--no_explain` | Flag | off (explain on) | Suppress the plain-language run summary. It is **on** by default: printed at the end and saved as `explanation.txt`. |
| `--no_print_info` | Flag | off (info on) | Skip the explanatory banner printed before a run. |

#### Configuration Tools

| Argument | Type | Default | Description |
| --- | --- | --- | --- |
| `--config_path` | String | `None` | JSON configuration file to read. |
| `--generate_config` | String | `None` | Write a template naming every key to this path, then exit. |
| `--config_preset` | String | `default` | Preset to inject into that template: `default`, `lima` or `arequipa`. |

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
* **Explanation:** The run in plain language, printed and saved. See §3's output list.