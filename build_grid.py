#!/usr/bin/env python3
"""
build_grid.py  -  real gridded paddock condition for several properties.

Per grid cell, indexed by the ACTUAL satellite pass date (newest first, up to 5):
  cover_dN / green_dN / coverdate_dN   DEA fractional cover for that pass
  temp_dN  / tempdate_dN               Landsat surface temperature for that pass
Plus weekly rainfall (SILO daily totalled per week), rain_wN with weeklabel_wN.

Each date slot is one real acquisition (tiles from the same day mosaicked across
the property). No carry-forward, so stepping the slider always moves to a genuinely
different satellite pass, labelled with the day it was taken. A cell is blank for a
date only if that pass did not image it.

cover = 100 - bare. green = photosynthetic fraction. dry = cover - green. In arid
country cover is mostly dry material, so the viewer also shows the green/dry/bare split.

Writes per property key into app_data/: <key>_grid.geojson, <key>_property.geojson
Runs on GitHub Actions. No invented numbers.
"""
import os, traceback
from collections import OrderedDict
from datetime import date, timedelta

import requests
import geopandas as gpd
from shapely.geometry import shape, box
from shapely.ops import unary_union
from pystac_client import Client
import planetary_computer as pc
import rasterio
from rasterio.crs import CRS
from rasterstats import zonal_stats

PROPERTIES = {
    "plevna":    {"name": "Plevna Downs",      "lots": ["6SP196201"]},
    "tilbooroo": {"name": "Tilbooroo Station", "lots": ["1336PH1989"]},
}
CELL_M   = int(os.environ.get("CELL_M", "1000"))
OUT      = "app_data"; os.makedirs(OUT, exist_ok=True)
ALBERS   = 3577
END      = date.today()
WEEKS    = [(END - timedelta(days=7*(i+1)), END - timedelta(days=7*i)) for i in range(4)]
LOOKBACK = 70                       # days back to gather satellite passes
MAX_DATES = 5                       # how many recent pass dates to keep per layer
MAX_SCENES = 16                     # safety cap on scene reads per layer
SILO_DAILY = "https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/daily/daily_rain/"
WGS84 = CRS.from_epsg(4326)

def log(*a): print(*a, flush=True)
def to_https(h):
    if h.startswith("s3://"):
        b, _, k = h[5:].partition("/"); return f"https://{b}.s3.ap-southeast-2.amazonaws.com/{k}"
    return h

def cadastre(lots):
    where = "lotplan IN (%s)" % ",".join(f"'{l}'" for l in lots)
    url = "https://spatial-gis.information.qld.gov.au/arcgis/rest/services/Basemaps/FoundationData/FeatureServer/2/query"
    r = requests.get(url, params=dict(where=where, outFields="lotplan",
                     returnGeometry="true", outSR=4326, f="geojson"), timeout=90)
    r.raise_for_status()
    feats = r.json().get("features", [])
    if not feats: raise RuntimeError("no parcels for " + str(lots))
    return unary_union([shape(f["geometry"]) for f in feats])

def make_grid(prop_al):
    minx, miny, maxx, maxy = prop_al.total_bounds
    cells = []; x = minx
    while x < maxx:
        y = miny
        while y < maxy:
            cells.append(box(x, y, x + CELL_M, y + CELL_M)); y += CELL_M
        x += CELL_M
    g = gpd.GeoDataFrame(geometry=cells, crs=ALBERS)
    whole = prop_al.union_all() if hasattr(prop_al, "union_all") else prop_al.unary_union
    g = g[g.intersects(whole)].copy()
    g["geometry"] = g.intersection(whole)
    g = g[~g.is_empty].reset_index(drop=True)
    g["cell_id"] = [f"C{i:03d}" for i in range(len(g))]
    return g

def zvals(grid, href, crs, nodata):
    gg = grid.to_crs(crs)
    st = zonal_stats(gg, f"/vsicurl/{to_https(href)}", stats=["mean"], nodata=nodata, all_touched=True)
    return [None if s["mean"] is None else s["mean"] for s in st]

dea = Client.open("https://explorer.dea.ga.gov.au/stac")
mpc = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)

def cover_by_date(grid, bbox, n):
    """OrderedDict date(newest first) -> (cover_list, green_list), tiles mosaicked."""
    items = sorted(dea.search(collections=["ga_ls_fc_3"], bbox=bbox,
                   datetime=f"{END-timedelta(days=LOOKBACK)}/{END}").items(),
                   key=lambda x: x.datetime, reverse=True)
    out = OrderedDict(); reads = 0
    for it in items:
        d = it.datetime.date()
        if d not in out and len(out) >= MAX_DATES: continue
        try:
            crs = rasterio.open(f"/vsicurl/{to_https(it.assets['bs'].href)}").crs
            bs = zvals(grid, it.assets["bs"].href, crs, 255); reads += 1
            if not any(v is not None for v in bs): continue
            pv = zvals(grid, it.assets["pv"].href, crs, 255)
            cov = [None if b is None else round(100-b) for b in bs]
            grn = [None if p is None else round(p) for p in pv]
            if d not in out: out[d] = (cov[:], grn[:])
            else:
                C, G = out[d]
                for j in range(n):
                    if C[j] is None and cov[j] is not None: C[j] = cov[j]; G[j] = grn[j]
            if reads >= MAX_SCENES: break
        except Exception as e:
            log("  cover scene skipped:", e)
    return out

def temp_by_date(grid, bbox, n):
    items = sorted(mpc.search(collections=["landsat-c2-l2"], bbox=bbox,
                   datetime=f"{END-timedelta(days=LOOKBACK)}/{END}",
                   query={"eo:cloud_cover": {"lt": 60}}).items(),
                   key=lambda x: x.datetime, reverse=True)
    out = OrderedDict(); reads = 0
    for it in items:
        d = it.datetime.date()
        if d not in out and len(out) >= MAX_DATES: continue
        try:
            crs = rasterio.open(f"/vsicurl/{it.assets['lwir11'].href}").crs
            dn = zvals(grid, it.assets["lwir11"].href, crs, 0); reads += 1
            if not any(v is not None for v in dn): continue
            t = [None if v is None else round(v*0.00341802 + 149.0 - 273.15, 1) for v in dn]
            if d not in out: out[d] = t
            else:
                T = out[d]
                for j in range(n):
                    if T[j] is None and t[j] is not None: T[j] = t[j]
            if reads >= MAX_SCENES: break
        except Exception as e:
            log("  temp scene skipped:", e)
    return out

for key, meta in PROPERTIES.items():
    try:
        log("===", key, meta["lots"])
        prop = gpd.GeoDataFrame(geometry=[cadastre(meta["lots"])], crs=4326)
        prop.to_file(f"{OUT}/{key}_property.geojson", driver="GeoJSON")
        prop_al = prop.to_crs(ALBERS)
        grid = make_grid(prop_al); n = len(grid)
        bbox = [float(v) for v in prop.total_bounds]
        log(f"{key}: {n} cells, {prop_al.area.iloc[0]/1e4:,.0f} ha")

        cov = cover_by_date(grid, bbox, n)
        for i, d in enumerate(list(cov.keys())[:MAX_DATES]):
            C, G = cov[d]
            grid[f"cover_d{i}"] = C; grid[f"green_d{i}"] = G; grid[f"coverdate_d{i}"] = str(d)
        log(f"{key}: cover dates {list(cov.keys())[:MAX_DATES]}")

        tmp = temp_by_date(grid, bbox, n)
        for i, d in enumerate(list(tmp.keys())[:MAX_DATES]):
            grid[f"temp_d{i}"] = tmp[d]; grid[f"tempdate_d{i}"] = str(d)
        log(f"{key}: temp dates {list(tmp.keys())[:MAX_DATES]}")

        for i, (ws, we) in enumerate(WEEKS):
            total = None; days = 0; d = ws
            while d < we:
                u = f"{SILO_DAILY}{d.year}/{d.year}{d.month:02d}{d.day:02d}.daily_rain.tif"
                try:
                    vals = zvals(grid, u, WGS84, None)
                    if total is None: total = [0.0]*n
                    for j, v in enumerate(vals):
                        if v is not None: total[j] += v
                    days += 1
                except Exception:
                    pass
                d += timedelta(days=1)
            grid[f"rain_w{i}"] = [None if total is None else round(total[j], 1) for j in range(n)]
            grid[f"weeklabel_w{i}"] = f"{ws:%d %b}-{we:%d %b}"
            log(f"{key} w{i} rain: {days} days")

        grid.to_crs(4326).to_file(f"{OUT}/{key}_grid.geojson", driver="GeoJSON")
        log(f"{key}: wrote grid")
    except Exception as e:
        log(key, "FAILED:", e); traceback.print_exc()

log("Done.")
