#!/usr/bin/env python
"""
Downloads road geometry from OpenStreetMap, for drawing on a search map.

Access is the question a site count cannot answer. A canyon wall that accepts a hundred
detectors is worth something quite different depending on whether a road runs along the
rim or the nearest track is forty kilometres away, and until now the maps this project
writes said nothing about it at all.

**This fetches geometry to draw, not a distance raster to screen on.** The searcher has
a separate, older facility for the latter -- ``road_map_path`` takes an aligned
distance-to-road GeoTIFF and ``max_road_dist_km`` cuts on it -- which nothing has ever
produced and no configuration uses. Drawing came first deliberately: it costs no
distance transform, it makes no claim about what is deployable, and it answers the
question a reader actually has when looking at a map.

    oroscope-fetch-roads --dem input/dem/arequipa_SRTMGL1.tif
    oroscope-fetch-roads --bbox -17.39 -14.56 -73.59 -70.09 --out input/roads/arequipa.geojson

Data is © OpenStreetMap contributors, ODbL. That attribution belongs on any figure that
shows it, and the loader carries it in the file so it cannot be separated from the data.

Queried through Overpass, which is a shared free service:

- The bbox is **tiled** and fetched a piece at a time. A single query over a 3.5-degree
  box times the server out with a 504; nine smaller ones succeed.
- There is a pause between requests, and one retry per tile against a second mirror.
  Neither is optional politeness -- Overpass rate-limits, and hammering it gets the
  address blocked.
- Only the classes worth drawing are asked for. Every residential street in Arequipa is
  a lot of geometry to render behind a mask and tells a reader nothing.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

# Ordered coarse to fine. A map at this scale wants the network that reaches a region,
# not every street inside a town.
ROAD_CLASSES = ("motorway", "trunk", "primary", "secondary", "tertiary")

# Mirrors, tried in order. The main instance is frequently overloaded and answers a
# perfectly good query with a 504.
OVERPASS_MIRRORS = (
    "https://overpass.kumi.systems/api/interpreter",
    "https://overpass-api.de/api/interpreter",
    "https://overpass.osm.ch/api/interpreter",
)

USER_AGENT = "oroscope/0.1 (terrain site search; https://github.com/mbustama/oroscope)"
ATTRIBUTION = "© OpenStreetMap contributors, ODbL (https://www.openstreetmap.org/copyright)"


def dem_bounds(dem_path):
    """
    South, north, west and east edges of a DEM, in degrees.

    Parameters
    ----------
    dem_path : str
        Path to the GeoTIFF.

    Returns
    -------
    tuple of float
        ``(south, north, west, east)``.
    """
    import tifffile as tiff
    from oroscope import site_searcher as ss

    cell_deg, rows, _ = ss.read_dem_geometry(dem_path)
    lat, lon = ss.read_dem_origin(dem_path)
    if cell_deg is None or lat is None:
        raise SystemExit(f"{dem_path} carries no pixel scale or tiepoint; use --bbox")
    with tiff.TiffFile(dem_path) as handle:
        rows, cols = handle.pages[0].shape[:2]
    return (lat - rows * cell_deg, lat, lon, lon + cols * cell_deg)


def tiles(bounds, step_deg=1.2):
    """
    Splits a bounding box into pieces small enough for one Overpass query.

    Parameters
    ----------
    bounds : tuple of float
        ``(south, north, west, east)`` in degrees.
    step_deg : float, optional
        Maximum span of a tile in either direction.

    Returns
    -------
    list of tuple
        ``(south, north, west, east)`` pieces covering the box.

    Examples
    --------
    >>> from oroscope import fetch_roads
    >>> len(fetch_roads.tiles((-17.4, -14.6, -73.6, -70.1), step_deg=1.5))
    6
    >>> fetch_roads.tiles((-1.0, 0.0, 0.0, 1.0), step_deg=2.0)
    [(-1.0, 0.0, 0.0, 1.0)]
    """
    south, north, west, east = bounds
    out = []
    lat = south
    while lat < north - 1e-9:
        lat_hi = min(lat + step_deg, north)
        lon = west
        while lon < east - 1e-9:
            lon_hi = min(lon + step_deg, east)
            out.append((lat, lat_hi, lon, lon_hi))
            lon = lon_hi
        lat = lat_hi
    return out


def query_for(bounds, classes=ROAD_CLASSES):
    """The Overpass QL for one tile."""
    south, north, west, east = bounds
    pattern = "|".join(classes)
    return (f'[out:json][timeout:120];'
            f'way["highway"~"^({pattern})$"]'
            f'({south:.5f},{west:.5f},{north:.5f},{east:.5f});'
            f'out geom;')


def fetch_tile(bounds, classes=ROAD_CLASSES, timeout=180, mirrors=OVERPASS_MIRRORS):
    """
    Fetches one tile, trying each mirror in turn.

    Parameters
    ----------
    bounds : tuple of float
        ``(south, north, west, east)``.
    classes : sequence of str, optional
        Highway classes to ask for.
    timeout : int, optional
        Seconds to wait for a response.
    mirrors : sequence of str, optional
        Overpass endpoints, tried in order.

    Returns
    -------
    list of dict
        Overpass ``way`` elements, each carrying a ``geometry`` list.

    Raises
    ------
    RuntimeError
        If every mirror failed.
    """
    payload = urllib.parse.urlencode({"data": query_for(bounds, classes)}).encode()
    problems = []
    for url in mirrors:
        request = urllib.request.Request(url, data=payload,
                                         headers={"User-Agent": USER_AGENT})
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return json.loads(response.read()).get("elements", [])
        except (urllib.error.URLError, urllib.error.HTTPError, ValueError, OSError) as e:
            problems.append(f"{url.split('/')[2]}: {type(e).__name__}")
            time.sleep(2.0)
    raise RuntimeError("every Overpass mirror failed -- " + "; ".join(problems))


PLACE_TYPES = ("city", "town", "village")


def places_query_for(bounds, types=PLACE_TYPES):
    """The Overpass QL for populated places in one tile."""
    south, north, west, east = bounds
    pattern = "|".join(types)
    return (f'[out:json][timeout:120];'
            f'node["place"~"^({pattern})$"]'
            f'({south:.5f},{west:.5f},{north:.5f},{east:.5f});'
            f'out body;')


def places_to_geojson(elements, types=PLACE_TYPES):
    """
    Turns Overpass place nodes into a GeoJSON FeatureCollection.

    Keeps the name, the place class and the population where OSM has one, which is what
    decides whether a place is worth a label on a crowded map.
    """
    features, seen = [], set()
    for element in elements:
        if element.get("id") in seen or "lat" not in element:
            continue
        seen.add(element.get("id"))
        tags = element.get("tags") or {}
        name = tags.get("name")
        if not name:
            continue
        try:
            population = int(str(tags.get("population", "")).replace(",", ""))
        except ValueError:
            population = None
        features.append({
            "type": "Feature",
            "properties": {"name": name, "place": tags.get("place"),
                           "population": population},
            "geometry": {"type": "Point",
                         "coordinates": [round(element["lon"], 5),
                                         round(element["lat"], 5)]},
        })
    return {"type": "FeatureCollection", "attribution": ATTRIBUTION,
            "types": list(types), "features": features}


def load_places(path):
    """
    Reads a places file written by this tool.

    Parameters
    ----------
    path : str
        Path to the GeoJSON.

    Returns
    -------
    list
        ``(latitude, longitude, name, place_class, population)`` tuples, most important
        first -- cities before towns before villages, and within a class the larger
        population first. A map that can only label a dozen places should label those.
        Empty when the file is absent or unreadable.
    """
    if not path or not os.path.exists(path):
        return []
    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, ValueError):                    # pragma: no cover - defensive
        return []
    rank = {"city": 0, "town": 1, "village": 2}
    out = []
    for feature in data.get("features", []):
        lon, lat = (feature.get("geometry") or {}).get("coordinates", (None, None))
        props = feature.get("properties") or {}
        if lat is None or not props.get("name"):
            continue
        out.append((lat, lon, props["name"], props.get("place"),
                    props.get("population")))
    out.sort(key=lambda p: (rank.get(p[3], 9), -(p[4] or 0)))
    return out


def fetch_places(bounds, out_path, types=PLACE_TYPES, step_deg=1.2, pause=1.5,
                 quiet=False):
    """
    Fetches populated places over a bounding box and writes one GeoJSON.

    Same tiling and pacing as :func:`fetch`, for the same reasons.
    """
    pieces = tiles(bounds, step_deg)
    elements = []
    for i, piece in enumerate(pieces, 1):
        if not quiet:
            print(f"   places tile {i}/{len(pieces)} ...", end=" ", flush=True)
        payload = urllib.parse.urlencode(
            {"data": places_query_for(piece, types)}).encode()
        got = []
        problems = []
        for url in OVERPASS_MIRRORS:
            request = urllib.request.Request(url, data=payload,
                                             headers={"User-Agent": USER_AGENT})
            try:
                with urllib.request.urlopen(request, timeout=180) as response:
                    got = json.loads(response.read()).get("elements", [])
                break
            except (urllib.error.URLError, urllib.error.HTTPError,
                    ValueError, OSError) as e:
                problems.append(f"{url.split('/')[2]}: {type(e).__name__}")
                time.sleep(2.0)
        else:                                        # pragma: no cover - all mirrors down
            raise RuntimeError("every Overpass mirror failed -- " + "; ".join(problems))
        elements += got
        if not quiet:
            print(f"{len(got)} places")
        if i < len(pieces):
            time.sleep(pause)

    collection = places_to_geojson(elements, types)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(collection, f, indent=1)
    if not quiet:
        print(f"\n   {len(collection['features']):,} places -> {out_path}")
    return collection


def to_geojson(elements, classes=ROAD_CLASSES):
    """
    Turns Overpass ways into a GeoJSON FeatureCollection, keeping only what is drawn.

    Overpass JSON is verbose and carries every tag; a map needs the geometry, the class
    and a name. Trimming here is what keeps a region's roads a few megabytes rather than
    a few tens.

    Parameters
    ----------
    elements : iterable of dict
        Overpass ``way`` elements.
    classes : sequence of str, optional
        Recorded in the collection's properties, so a file says what was asked for --
        an empty result then means "no roads of these classes here" rather than
        "something went wrong".

    Returns
    -------
    dict
        A GeoJSON FeatureCollection of LineStrings, with ``highway`` and ``name``.
    """
    features, seen = [], set()
    for element in elements:
        geometry = element.get("geometry") or []
        if len(geometry) < 2:
            continue
        if element.get("id") in seen:               # tiles overlap at their seams
            continue
        seen.add(element.get("id"))
        tags = element.get("tags") or {}
        features.append({
            "type": "Feature",
            "properties": {"highway": tags.get("highway"),
                           "name": tags.get("name") or tags.get("ref")},
            "geometry": {"type": "LineString",
                         "coordinates": [[round(p["lon"], 5), round(p["lat"], 5)]
                                         for p in geometry]},
        })
    return {"type": "FeatureCollection",
            "attribution": ATTRIBUTION,
            "classes": list(classes),
            "features": features}


def load_roads(path):
    """
    Reads a road file written by this tool.

    Parameters
    ----------
    path : str
        Path to the GeoJSON.

    Returns
    -------
    dict or None
        ``{"attribution", "classes", "lines"}`` where ``lines`` is a list of
        ``(highway_class, [(latitude, longitude), ...])``. ``None`` when the file is
        absent or unreadable, so a missing road file leaves a map without roads rather
        than without a map.
    """
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            data = json.load(f)
    except (OSError, ValueError):                    # pragma: no cover - defensive
        return None
    lines = []
    for feature in data.get("features", []):
        coords = (feature.get("geometry") or {}).get("coordinates") or []
        if len(coords) < 2:
            continue
        klass = (feature.get("properties") or {}).get("highway")
        lines.append((klass, [(lat, lon) for lon, lat in coords]))
    return {"attribution": data.get("attribution", ATTRIBUTION),
            "classes": data.get("classes", list(ROAD_CLASSES)),
            "lines": lines}


def fetch(bounds, out_path, classes=ROAD_CLASSES, step_deg=1.2, pause=1.5, quiet=False):
    """
    Fetches every tile of a bounding box and writes one GeoJSON.

    Parameters
    ----------
    bounds : tuple of float
        ``(south, north, west, east)`` in degrees.
    out_path : str
        Where to write the GeoJSON.
    classes : sequence of str, optional
        Highway classes to ask for.
    step_deg : float, optional
        Tile size. Larger is fewer requests and more 504s.
    pause : float, optional
        Seconds between requests. Overpass is a shared free service.
    quiet : bool, optional
        Suppress progress.

    Returns
    -------
    dict
        The GeoJSON written.
    """
    pieces = tiles(bounds, step_deg)
    elements = []
    for i, piece in enumerate(pieces, 1):
        if not quiet:
            print(f"   tile {i}/{len(pieces)} "
                  f"({piece[0]:.2f}..{piece[1]:.2f}, {piece[2]:.2f}..{piece[3]:.2f}) ...",
                  end=" ", flush=True)
        got = fetch_tile(piece, classes)
        elements += got
        if not quiet:
            print(f"{len(got)} ways")
        if i < len(pieces):
            time.sleep(pause)

    collection = to_geojson(elements, classes)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(collection, f)
    if not quiet:
        size = os.path.getsize(out_path) / 1024.0
        print(f"\n   {len(collection['features']):,} roads -> {out_path} ({size:,.0f} KiB)")
        print(f"   {ATTRIBUTION}")
    return collection


def main():
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dem", type=str, default=None,
                        help="Take the bounding box from this DEM, which is what you "
                             "want when the roads are to be drawn over its search.")
    parser.add_argument("--bbox", type=float, nargs=4, default=None,
                        metavar=("SOUTH", "NORTH", "WEST", "EAST"),
                        help="Bounding box in degrees, instead of --dem.")
    parser.add_argument("--out", type=str, default=None,
                        help="Output GeoJSON (default: input/roads/<dem name>.geojson).")
    parser.add_argument("--classes", type=str, nargs="+", default=list(ROAD_CLASSES),
                        help=f"Highway classes to fetch (default: {' '.join(ROAD_CLASSES)}).")
    parser.add_argument("--places", action="store_true",
                        help="Also fetch populated places (city, town, village) into a "
                             "sibling _places.geojson, for marking towns on the map.")
    parser.add_argument("--places_only", action="store_true",
                        help="Fetch only the places, skipping the roads.")
    parser.add_argument("--step_deg", type=float, default=1.2,
                        help="Tile size in degrees. Larger is fewer requests, and more "
                             "gateway timeouts.")
    args = parser.parse_args()

    if args.dem:
        if not os.path.exists(args.dem):
            raise SystemExit(f"no such DEM: {args.dem}")
        bounds = dem_bounds(args.dem)
        default_name = os.path.splitext(os.path.basename(args.dem))[0] + ".geojson"
    elif args.bbox:
        bounds = tuple(args.bbox)
        default_name = "roads.geojson"
    else:
        raise SystemExit("give --dem or --bbox")

    out = args.out or os.path.join("input", "roads", default_name)
    print(f"bbox: {bounds[0]:.4f}..{bounds[1]:.4f} lat, {bounds[2]:.4f}..{bounds[3]:.4f} lon")
    if not args.places_only:
        print(f"classes: {', '.join(args.classes)}\n")
        fetch(bounds, out, args.classes, args.step_deg)
    if args.places or args.places_only:
        places_out = out.replace(".geojson", "_places.geojson")
        print(f"\nplaces: {', '.join(PLACE_TYPES)}\n")
        fetch_places(bounds, places_out, step_deg=args.step_deg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
