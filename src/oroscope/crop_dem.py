#!/usr/bin/env python
"""
Cuts a lat/lon window out of a DEM into a smaller GeoTIFF.

Two reasons this exists. A regional study wants one canyon or one plateau rather than
the whole 10204x12603 tile the download gives you. And comparing experiments requires
them to cover *identical* ground: the combiner overlays masks pixel for pixel, so
GRAND and TAMBO have to be run against the same crop, not two crops that happen to
overlap.

The window is snapped outward to whole pixels and the tiepoint recomputed for the
crop's own corner, so the output is georeferenced in its own right.

    python crop_dem.py ../input/dem/arequipa_SRTMGL1.tif ../input/dem/colca.tif \
        --north -15.30 --south -15.80 --west -72.40 --east -71.60
"""

from __future__ import annotations

import argparse
import os

import numpy as np

__all__ = ["crop", "read_geo", "main"]
import tifffile as tiff


def read_geo(path: str) -> tuple[float, float, float, float, int, int]:
    """
    Pixel size in degrees and the north-west corner, from the GeoTIFF tags.

    Parameters
    ----------
    path : str
        Path to a geographic GeoTIFF carrying ``ModelPixelScaleTag`` and
        ``ModelTiepointTag``.

    Returns
    -------
    tuple
        ``(cell_x_deg, cell_y_deg, lon0, lat0, rows, cols)``, with ``lon0``/``lat0``
        the north-west corner.
    """
    with tiff.TiffFile(path) as tf:
        page = tf.pages[0]
        scale = page.tags["ModelPixelScaleTag"].value
        tie = page.tags["ModelTiepointTag"].value
        rows, cols = int(page.shape[0]), int(page.shape[1])
    return float(scale[0]), float(scale[1]), float(tie[3]), float(tie[4]), rows, cols


def crop(src: str, dst: str, north: float, south: float, west: float,
         east: float) -> dict:
    """
    Cuts a geographic window out of a DEM, writing a GeoTIFF that stands on its own.

    The point of cropping is not disk space, it is sampling: a department at
    ``downsample_factor`` 4 and ``candidate_stride`` 5 costs area and fragments a thin
    mask, while a crop small enough to run at 1 and 1 does neither. Every unbiased
    number this project quotes comes from a crop made here.

    The window is snapped outward to whole pixels, so the result contains the requested
    box rather than approximating it, and is clipped to the DEM. A window that misses
    the DEM entirely raises rather than writing an empty file.

    **The crop carries its own tiepoint**, not the parent's. That is what lets it be
    searched with no reference back to where it came from — and it is why a crop's
    ``origin_lat``/``origin_lon`` differ from the parent's, which is correct and not
    drift.

    Parameters
    ----------
    src : str
        The GeoTIFF to cut from.
    dst : str
        Where to write the crop.
    north, south : float
        Latitude bounds in degrees. ``north`` is the larger (less negative) value.
    west, east : float
        Longitude bounds in degrees.

    Returns
    -------
    dict
        ``path``, the grid as ``rows`` and ``cols``, the crop's own ``origin_lat`` and
        ``origin_lon`` with the matching ``south`` and ``east`` edges,
        ``cell_size_deg``, and the elevation range as ``z_min`` and ``z_max``.

    Raises
    ------
    SystemExit
        If the requested window does not overlap the DEM, reporting what the DEM
        actually covers.

    Examples
    --------
    Cutting the Colca crop out of the Arequipa department DEM, which is how
    ``input/dem/colca.tif`` was made::

        from oroscope import crop_dem

        info = crop_dem.crop("input/dem/arequipa_SRTMGL1.tif", "input/dem/colca.tif",
                             north=-15.30, south=-15.85, west=-72.40, east=-71.55)
        print(info["rows"], info["cols"], info["z_min"], info["z_max"])
    """
    cell_x_deg, cell_y_deg, lon0, lat0, rows, cols = read_geo(src)

    # Rows run north to south, columns west to east
    r_start = int(np.floor((lat0 - north) / cell_y_deg))
    r_stop = int(np.ceil((lat0 - south) / cell_y_deg))
    c_start = int(np.floor((west - lon0) / cell_x_deg))
    c_stop = int(np.ceil((east - lon0) / cell_x_deg))

    r_start, c_start = max(0, r_start), max(0, c_start)
    r_stop, c_stop = min(rows, r_stop), min(cols, c_stop)
    if r_stop <= r_start or c_stop <= c_start:
        raise SystemExit(
            f"requested window does not overlap the DEM.\n"
            f"  DEM covers lat {lat0:.4f} .. {lat0 - rows * cell_y_deg:.4f}, "
            f"lon {lon0:.4f} .. {lon0 + cols * cell_x_deg:.4f}")

    with tiff.TiffFile(src) as tf:
        data = tf.pages[0].asarray()
    window = np.ascontiguousarray(data[r_start:r_stop, c_start:c_stop])

    # The crop's own north-west corner, so it stands alone as a georeferenced file
    out_lat = lat0 - r_start * cell_y_deg
    out_lon = lon0 + c_start * cell_x_deg
    tiff.imwrite(
        dst, window,
        extratags=[
            (33550, "d", 3, (cell_x_deg, cell_y_deg, 0.0)),
            (33922, "d", 6, (0.0, 0.0, 0.0, out_lon, out_lat, 0.0)),
        ],
    )
    return dict(path=dst, rows=window.shape[0], cols=window.shape[1],
                origin_lat=out_lat, origin_lon=out_lon,
                south=out_lat - window.shape[0] * cell_y_deg,
                east=out_lon + window.shape[1] * cell_x_deg,
                cell_size_deg=cell_y_deg,
                z_min=float(np.nanmin(window)), z_max=float(np.nanmax(window)))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--north", type=float, required=True)
    ap.add_argument("--south", type=float, required=True)
    ap.add_argument("--west", type=float, required=True)
    ap.add_argument("--east", type=float, required=True)
    args = ap.parse_args()

    if args.north <= args.south:
        raise SystemExit("--north must be greater than --south")
    if args.east <= args.west:
        raise SystemExit("--east must be greater than --west")

    info = crop(args.src, args.dst, args.north, args.south, args.west, args.east)
    size_mb = os.path.getsize(info["path"]) / 1e6
    print(f"wrote {info['path']}  ({size_mb:.1f} MB)")
    print(f"   {info['rows']} x {info['cols']} px at {info['cell_size_deg']*3600:.2f} arcsec")
    print(f"   lat {info['origin_lat']:.4f} .. {info['south']:.4f}")
    print(f"   lon {info['origin_lon']:.4f} .. {info['east']:.4f}")
    print(f"   elevation {info['z_min']:.0f} .. {info['z_max']:.0f} m")
    print(f"\n   origin_lat / origin_lon for the config: "
          f"{info['origin_lat']:.6f}, {info['origin_lon']:.6f}")


if __name__ == "__main__":
    main()
