#!/usr/bin/env python3
"""
build_grid.py  -  real gridded paddock condition, weekly, for several properties.

For each property and each of the last 4 weekly windows (w0 = this week back to
w3 = 3 weeks ago), per grid cell:
  cover_wN / green_wN / coverdate_wN   DEA fractional cover
  temp_wN  / tempdate_wN               Landsat surface temperature
  rain_wN                              SILO daily rainfall totalled over the week

Cover and temperature use "latest clear scene on or before each week" (last
observation carried forward), so a week with no satellite pass still shows the
most recent reading, dated, instead of a blank. Scenes that do not actually
overlap the property (a neighbouring Landsat tile) are rejected.

Writes per property key into app_data/: <key>_grid.geojson, <key>_property.geojson
Runs on GitHub Actions. No invented numbers.
"""
import os, traceback
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
CELL_M  = int(os.environ.get("CELL_M", "1000"))
OUT     = "app_data"; os.makedirs(OUT, exist_ok=True)
ALBERS  = 3577
END     = date.today()
WEEKS   = [(END - timedelta(days=7*(i+1)), END - timedelta(days=7*i)) for i in range(4)]
LOOKBACK = 42                       # days back to search for a usable cover/temp scene
MIN_OVERLAP = 0.30                  # a scene must cover at least this fraction of cells
MAX_SCENES = 8                      # cap reads per layer per property
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
    poly = g.union_all() if hasattr(g, "union_all") else g.unary_union
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

def frac_valid(vals):
    return sum(1 for v in vals if v is not None) / max(1, len(vals))

dea = Client.open("https://explorer.dea.ga.gov.au/stac")
mpc = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)

def usable_cover(grid, bbox):
    """Newest-first list of (date, cover_list, green_list) for overlapping FC scenes."""
    items = sorted(dea.search(collections=["ga_ls_fc_3"], bbox=bbox,
                   datetime=f"{END-timedelta(days=LOOKBACK)}/{END}").items(),
                   key=lambda x: x.datetime, reverse=True)
    out = []
    for it in items:
        try:
            crs = rasterio.open(f"/vsicurl/{to_https(it.assets['bs'].href)}").crs
            bs = zvals(grid, it.assets["bs"].href, crs, 255)
            if frac_valid(bs) < MIN_OVERLAP:
                continue
            pv = zvals(grid, it.assets["pv"].href, crs, 255)
            out.append((it.datetime.date(),
                        [None if b is None else round(100 - b) for b in bs],
                        [None if p is None else round(p) for p in pv]))
            if len(out) >= MAX_SCENES: break
        except Exception as e:
            log("  cover scene skipped:", e)
    return out

def usable_temp(grid, bbox):
    """Newest-first list of (date, temp_list) for overlapping, mostly-clear Landsat."""
    items = sorted(mpc.search(collections=["landsat-c2-l2"], bbox=bbox,
                   datetime=f"{END-timedelta(days=LOOKBACK)}/{END}",
                   query={"eo:cloud_cover": {"lt": 60}}).items(),
                   key=lambda x: x.datetime, reverse=True)
    out = []
    for it in items:
        try:
            crs = rasterio.open(f"/vsicurl/{it.assets['lwir11'].href}").crs
            dn = zvals(grid, it.assets["lwir11"].href, crs, 0)
            if frac_valid(dn) < MIN_OVERLAP:
                continue
            out.append((it.datetime.date(),
                        [None if d is None else round(d*0.00341802 + 149.0 - 273.15, 1) for d in dn]))
            if len(out) >= MAX_SCENES: break
        except Exception as e:
            log("  temp scene skipped:", e)
    return out

for key, meta in PROPERTIES.items():
    try:
        log("===", key, meta["lots"])
        prop = gpd.GeoDataFrame(geometry=[cadastre(meta["lots"])], crs=4326)
        prop.to_file(f"{OUT}/{key}_property.geojson", driver="GeoJSON")
        prop_al = prop.to_crs(ALBERS)
        grid = make_grid(prop_al)
        bbox = [float(v) for v in prop.total_bounds]
        log(f"{key}: {len(grid)} cells, {prop_al.area.iloc[0]/1e4:,.0f} ha")

        cov = usable_cover(grid, bbox)
        tmp = usable_temp(grid, bbox)
        log(f"{key}: {len(cov)} cover scenes, {len(tmp)} temp scenes available")

        for i, (ws, we) in enumerate(WEEKS):
            tag = f"w{i}"
            # cover: latest usable scene on or before this week's end
            pick = next((c for c in cov if c[0] <= we), None)
            if pick:
                grid[f"cover_{tag}"] = pick[1]; grid[f"green_{tag}"] = pick[2]; grid[f"coverdate_{tag}"] = str(pick[0])
            else:
                grid[f"cover_{tag}"] = None; grid[f"green_{tag}"] = None; grid[f"coverdate_{tag}"] = None
            # temp: latest usable scene on or before this week's end
            pick = next((t for t in tmp if t[0] <= we), None)
            if pick:
                grid[f"temp_{tag}"] = pick[1]; grid[f"tempdate_{tag}"] = str(pick[0])
            else:
                grid[f"temp_{tag}"] = None; grid[f"tempdate_{tag}"] = None
            # rain: SILO daily totalled strictly within this week
            try:
                total = None; days = 0; d = ws
                while d < we:
                    u = f"{SILO_DAILY}{d.year}/{d.year}{d.month:02d}{d.day:02d}.daily_rain.tif"
                    try:
                        vals = zvals(grid, u, WGS84, None)
                        if total is None: total = [0.0]*len(vals)
                        for j, v in enumerate(vals):
                            if v is not None: total[j] += v
                        days += 1
                    except Exception:
                        pass
                    d += timedelta(days=1)
                grid[f"rain_{tag}"] = [None if total is None else round(total[j], 1) for j in range(len(grid))]
                log(f"{key} {tag} rain: {days} days")
            except Exception as e:
                log(f"{key} {tag} rain failed:", e); grid[f"rain_{tag}"] = None

        grid.to_crs(4326).to_file(f"{OUT}/{key}_grid.geojson", driver="GeoJSON")
        log(f"{key}: wrote grid")
    except Exception as e:
        log(key, "FAILED:", e); traceback.print_exc()

log("Done.")
