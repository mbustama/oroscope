
"""
Downloads the elevation models a search needs, and writes configurations for them.

Fetches the bundled regions from OpenTopography into ``input/dem/`` and generates a
ready-to-run JSON config for each in ``config/``::

    oroscope-fetch-dem --region arequipa --open_topography_api_key KEY

Four regions are defined. Three departments at 1 arc-second from **the same dataset**,
so that runs over them are comparable -- ``arequipa`` (129 Mpx), ``lima`` (105 Mpx) and
``ancash`` (69 Mpx), all SRTMGL1 -- and ``peru``, the whole country, at 3 arc-seconds
(SRTMGL3) because that is the only resolution that fits either a desktop's memory or the
API's own area limit. Omit ``--region`` to fetch all of them.

**Getting a key.** It is free and takes a minute. Register at
https://portal.opentopography.org/myopentopo, sign in, then open *myOpenTopo
Authorizations and API Key* from the account menu and copy the key. Pass it as
``--open_topography_api_key``, or set ``OPENTOPOGRAPHY_API_KEY`` in the environment,
which keeps it out of your shell history and out of any file that might be committed.

**Requests are capped by area,** and the cap is per dataset: 450,000 km² for every 30 m
dataset, 4,050,000 km² for the 90 m ones. That is why ``peru`` is SRTMGL3 -- its
bounding box is about 2.86 million km², six times over the 30 m limit.

For any other region, download the tiles from the OpenTopography portal, merge them
into one GeoTIFF if the area spans several, and cut the window you want with
:mod:`crop_dem`.

This was ``setup.py``, whose name made ``pip install`` run the downloader instead of
building the package.
"""

from __future__ import annotations
import argparse
import sys
import os
import urllib.request
import urllib.error
import json
from tqdm import tqdm

# ==========================================
#          UI THEME & FORMATTING
# ==========================================
if sys.platform == 'win32':
    os.system('') 

def supports_color():
    supported_platform = sys.platform != 'win32' or 'ANSICON' in os.environ or 'WT_SESSION' in os.environ or os.environ.get('TERM') == 'xterm-256color'
    is_a_tty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()
    return supported_platform and is_a_tty

def supports_emoji():
    if sys.stdout.encoding:
        return sys.stdout.encoding.lower() == 'utf-8'
    return False

USE_COLOR = supports_color()
USE_EMOJI = supports_emoji()

class C:
    HEADER = '\033[96m' if USE_COLOR else ''
    OK = '\033[92m' if USE_COLOR else ''
    WARN = '\033[93m' if USE_COLOR else ''
    FAIL = '\033[91m' if USE_COLOR else ''
    BOLD = '\033[1m' if USE_COLOR else ''
    MAGENTA = '\033[95m' if USE_COLOR else ''
    RESET = '\033[0m' if USE_COLOR else ''

class Icon:
    MAP = '🗺️  ' if USE_EMOJI else '[*] '
    GEAR = '⚙️  ' if USE_EMOJI else '[~] '
    DISK = '💾 ' if USE_EMOJI else '[S] '
    WARN = '⚠️  ' if USE_EMOJI else '[!] '
    CHECK = '✅ ' if USE_EMOJI else '[✓] '
    CROSS = '❌ ' if USE_EMOJI else '[x] '
    INFO = 'ℹ️  ' if USE_EMOJI else '[i] '

# ==========================================
#             TARGET REGIONS
# ==========================================
#: The regions ``oroscope-fetch-dem --region`` knows, each mapping to its dataset and
#: its bounding box in degrees. A ``#:`` comment rather than a plain one so that autodoc
#: picks it up: :doc:`data` cross-references this name, and module-level data with only
#: an ordinary comment above it is documented nowhere, leaving the reference to render
#: as unlinked text.
#:
#: Bounds are given as (West, East, South, North).
REGIONS = {
    # SRTMGL1, not the AW3D30 this used to fetch. The department runs are compared
    # against one another and a dataset difference would sit inside every comparison as
    # a confound -- Arequipa and Ancash are both SRTMGL1, so Lima is too. 98,222 km2 is
    # well inside the API's 450,000 km2 limit for a 30 m dataset. The older
    # lima_AW3D30.tif is left on disk for anything that referenced it.
    "lima": {
        "west": -78.07665824890137,
        "east": -75.39955615997313,
        "south": -13.252477566131276,
        "north": -10.228479499469358,
        "filename": "lima_SRTMGL1.tif",
        "preset": "lima",
        "demtype": "SRTMGL1"
    },
    "arequipa": {
        "west": -73.58612537384033,
        "east": -70.0852632522583,
        "south": -17.38995824658555,
        "north": -14.555380967667489,
        "filename": "arequipa_SRTMGL1.tif",
        "preset": "arequipa",
        "demtype": "SRTMGL1"
    },
    # Ancash, at 1 arc-second like Arequipa so the two are directly comparable. The
    # Cordillera Blanca and the Callejon de Huaylas: far steeper ground than Arequipa,
    # which is the point of running it. Bounds are OpenStreetMap's administrative
    # boundary for the department, queried from Nominatim rather than eyeballed.
    # 9,855 x 6,958 = 68.6 Mpx, 53% of Arequipa's 128.6, and 64,684 km2 -- comfortably
    # inside the 450,000 km2 the API allows for a 30 m dataset.
    "ancash": {
        "west": -78.6584805,
        "east": -76.7257441,
        "south": -10.7873076,
        "north": -8.0497090,
        "filename": "ancash_SRTMGL1.tif",
        "preset": "default",
        "demtype": "SRTMGL1"
    },
    # The whole country, at 3 arc-seconds rather than 1. Not a preference, and it is
    # forced twice over. By memory: 18.4 by 12.8 degrees is 3,052 Mpx at 1 arc-second
    # and 339 Mpx at 3, and only the latter fits in a desktop in one run -- see
    # config/grand_peru_survey.json for the arithmetic. And by the API, which caps a
    # request at "4,050,000 km2 for SRTM GL3, COP90 (90m resolution), and 450,000 km2
    # for all 30m resolution datasets" (opentopography.org/developers). This box is
    # ~2.86 million km2: inside the 90 m limit, six times over the 30 m one.
    #
    # A 90 m pixel is a survey instrument. It resolves the plateaux and ranges GRAND
    # deploys on; it does not resolve a canyon wall, so it is not the grid to run TAMBO
    # on. If SRTMGL3 is ever rejected, COP90 is the same resolution and the same limit.
    "peru": {
        "west": -81.4,
        "east": -68.6,
        "south": -18.4,
        "north": 0.0,
        "filename": "peru_SRTMGL3.tif",
        "preset": "default",
        "demtype": "SRTMGL3"
    }
}

# ==========================================
#               CORE LOGIC
# ==========================================
class TqdmUpTo(tqdm):
    """Callback class for urlretrieve to print a dynamic tqdm progress bar."""
    def update_to(self, b=1, bsize=1, tsize=None):
        if tsize is not None:
            self.total = tsize
        self.update(b * bsize - self.n)

def download_dem(region_name, bounds, api_key, output_dir):
    """
    Fetches one region's GeoTIFF from OpenTopography.

    Parameters
    ----------
    region_name : str
        A key of :data:`REGIONS`, naming the box and dataset to fetch.
    bounds : dict
        That region's entry: ``demtype`` and the ``south``/``north``/``west``/``east``
        bounds in degrees, plus the ``filename`` to write.
    api_key : str
        An OpenTopography key. Free, from
        https://portal.opentopography.org/myopentopo.
    output_dir : str
        Directory to write the ``.tif`` into. Created if absent.

    Returns
    -------
    str or None
        The path written, or ``None`` if the request failed.
    """
    url = (f"https://portal.opentopography.org/API/globaldem"
           f"?demtype={bounds['demtype']}"
           f"&south={bounds['south']}"
           f"&north={bounds['north']}"
           f"&west={bounds['west']}"
           f"&east={bounds['east']}"
           f"&outputFormat=GTiff"
           f"&API_Key={api_key}")
    
    filepath = os.path.join(output_dir, bounds['filename'])
    
    print(f"\n{C.BOLD}Processing Region: {region_name.title()}{C.RESET}")
    print(f"   {Icon.MAP}Requesting {bounds['demtype']} DEM from OpenTopography...")
    
    try:
        with TqdmUpTo(unit='B', unit_scale=True, unit_divisor=1024, miniters=1,
                      desc=f"      {Icon.INFO}Downloading", colour='magenta' if USE_COLOR else None) as t:
            urllib.request.urlretrieve(url, filepath, reporthook=t.update_to)
            
        print(f"      {Icon.CHECK}{C.OK}Saved to: {filepath}{C.RESET}")
        return filepath
    except urllib.error.HTTPError as e:
        print(f"      {Icon.CROSS}{C.FAIL}HTTP Error {e.code}: {e.reason}{C.RESET}")
        if e.code == 401:
            print(f"      {C.WARN}Make sure your OpenTopography API key is valid and active.{C.RESET}")
        return None
    except Exception as e:
        print(f"      {Icon.CROSS}{C.FAIL}Download failed: {e}{C.RESET}")
        return None

def generate_and_patch_config(region_name, preset, dem_filepath, config_dir=None):
    """
    Writes a default config for the region and points it at the DEM just downloaded.

    This used to shell out to ``site_searcher.py`` in the current directory, which was
    right when the modules were flat and is wrong now that they are a package: an
    installed copy has no such file anywhere, so the config step failed for every user
    who was not standing in ``src/``. :func:`site_searcher.generate_config` is the same
    code path the ``--generate_config`` flag takes.

    Parameters
    ----------
    region_name : str
        Key in :data:`REGIONS`; names the config file.
    preset : str
        Config preset to seed from, one of ``site_searcher.CONFIG_PRESETS``.
    dem_filepath : str
        The DEM to point ``dem_path`` at.
    config_dir : str, optional
        Where to write. Defaults to ``../config`` for continuity with the old layout.
    """
    config_dir = config_dir or os.path.join("..", "config")
    os.makedirs(config_dir, exist_ok=True)
    config_filename = os.path.join(config_dir, f"{region_name}_config.json")

    print(f"   {Icon.GEAR}Generating config file...")

    # 1. Generate the raw config through the library
    try:
        from oroscope import site_searcher as ss
        ss.generate_config(config_filename, preset)
    except Exception as e:
        print(f"      {Icon.CROSS}{C.FAIL}Failed to generate config: {e}{C.RESET}")
        return

    # 2. Patch the JSON to point to the newly downloaded DEM
    try:
        with open(config_filename, 'r') as f:
            config_data = json.load(f)
            
        # Use a forward-slash relative path starting from the config's directory
        rel_dem_path = os.path.relpath(dem_filepath, start=config_dir).replace('\\', '/')
        config_data['dem_path'] = rel_dem_path
        
        with open(config_filename, 'w') as f:
            json.dump(config_data, f, indent=4)
            
        print(f"      {Icon.CHECK}{C.OK}Config configured and saved to: {config_filename}{C.RESET}")
    except Exception as e:
        print(f"      {Icon.CROSS}{C.FAIL}Error patching JSON config: {e}{C.RESET}")

# ==========================================
#                  MAIN
# ==========================================
def main():
    """
    Command-line entry point.

    Downloads the DEM tiles a search needs and writes matching config files.

    Kept as a function so the console script declared in pyproject.toml has something
    to call.
    """
    parser = argparse.ArgumentParser(description="Oroscope - DEM download and configuration setup")
    parser.add_argument("--open_topography_api_key", type=str, default=None,
                        help="Your OpenTopography API key (Required to download DEMs). Free, from https://portal.opentopography.org/myopentopo -- register, then 'myOpenTopo Authorizations and API Key'. Read from the OPENTOPOGRAPHY_API_KEY environment variable if this flag is not given, which keeps it out of your shell history.")
    parser.add_argument("--region", type=str, default=None, choices=sorted(REGIONS),
                        help="Download one region rather than all of them. Peru is 302 MB and the whole set is over 550; asking for everything to get one is the wrong default and used to be the only option.")
    parser.add_argument("--output_dir", type=str, default=os.path.join("..", "input", "dem"),
                        help="Where the DEMs land. Default ../input/dem, which is right when standing in src/.")
    parser.add_argument("--config_dir", type=str, default=os.path.join("..", "config"),
                        help="Where the generated configs land. Default ../config.")

    args = parser.parse_args()
    api_key = args.open_topography_api_key or os.environ.get("OPENTOPOGRAPHY_API_KEY")

    print(f"\n{C.HEADER}===================================================={C.RESET}")
    print(f"{C.BOLD}   OROSCOPE OPENTOPOGRAPHY DEM SETUP{C.RESET}")
    print(f"{C.HEADER}===================================================={C.RESET}")

    # Explicit enforcement of the API Key
    if not api_key:
        print(f"\n{C.FAIL}{Icon.CROSS}ERROR: Missing Required Parameter.{C.RESET}")
        print(f"{C.WARN}Pass {C.BOLD}--open_topography_api_key{C.RESET}{C.WARN}, or set OPENTOPOGRAPHY_API_KEY in the environment.{C.RESET}")
        print("The key is free. Register at the link below, then open 'myOpenTopo Authorizations")
        print("and API Key' and copy the key from there. Without one you must download the TIF")
        print("files by hand from the OpenTopography portal.")
        print(f"Register for a free key at: {C.MAGENTA}https://portal.opentopography.org/myopentopo{C.RESET}\n")
        sys.exit(1)

    # Dependency check using generate_env.py. Against the installed module rather than a
    # site_searcher.py in the current directory: there used to be a hard exit here if
    # that file was not beside you, which made the whole tool unusable from anywhere but
    # the old flat src/ layout.
    try:
        from oroscope import generate_env, site_searcher
        print(f"\n{C.HEADER}===================================================={C.RESET}")
        print(f"{C.BOLD}   DEPENDENCY CHECK{C.RESET}")
        print(f"{C.HEADER}===================================================={C.RESET}")
        print(f"   {Icon.GEAR}Scanning the searcher for dependencies...")
        deps = generate_env.extract_dependencies(site_searcher.__file__)
        satisfied, missing = generate_env.check_installed_modules(deps)
        
        if missing:
            print(f"   {C.FAIL}{Icon.CROSS}Missing dependencies detected: {', '.join(missing)}{C.RESET}")
            print(f"   {C.WARN}It is highly recommended to run 'python generate_env.py' to setup your environment.{C.RESET}")
        else:
            print(f"   {C.OK}{Icon.CHECK}All script dependencies are satisfied!{C.RESET}")
    except ImportError:
        print(f"\n   {C.WARN}{Icon.WARN}Could not import 'generate_env.py' to verify dependencies. Skipping check.{C.RESET}")
    except Exception as e:
        print(f"\n   {C.FAIL}{Icon.CROSS}Error checking dependencies: {e}{C.RESET}")

    # Establish the ../input/dem/ directory structure
    print(f"\n{C.HEADER}===================================================={C.RESET}")
    print(f"{C.BOLD}   DOWNLOADING ASSETS{C.RESET}")
    print(f"{C.HEADER}===================================================={C.RESET}")
    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    print(f"   {Icon.DISK}Target directory verified: {C.MAGENTA}{os.path.abspath(output_dir)}{C.RESET}")

    wanted = [args.region] if args.region else list(REGIONS)
    written = []
    for region in wanted:
        bounds = REGIONS[region]
        downloaded_path = download_dem(region, bounds, api_key, output_dir)

        if downloaded_path:
            generate_and_patch_config(region, bounds['preset'], downloaded_path,
                                      args.config_dir)
            written.append(region)

    print(f"\n{C.HEADER}===================================================={C.RESET}")
    print(f"{C.OK}{Icon.CHECK}{C.BOLD}Setup Complete.{C.RESET}")
    if written:
        print("You can now run a search against a generated config:")
        for region in written:
            print(f"  {C.MAGENTA}oroscope --config_path "
                  f"{os.path.join(args.config_dir, region + '_config.json')}{C.RESET}")
        print("\nThe generated config is a default template. The tuned ones that produced")
        print("this project's published numbers live beside it in config/ -- for Peru,")
        print(f"  {C.MAGENTA}oroscope --config_path "
              f"{os.path.join(args.config_dir, 'grand_peru_survey.json')}{C.RESET}\n")


if __name__ == "__main__":
    main()
