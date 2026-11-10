#!/usr/bin/env python3
"""
build_grid.py  -  real gridded paddock condition, weekly, for several properties.

For each property and each of the last 4 weekly windows (w0 = this week back to
w3 = 3 weeks ago), computes per grid cell:
  cover_wN / green_wN / coverdate_wN   DEA fractional cover (latest scene in week)
  temp_wN  / tempdate_wN               Landsat surface temperature (latest in week)
  rain_wN                              SILO daily rainfall totalled over the week

Writes, per property key, into app_data/:
  <key>_grid.geojson      the grid with all weekly columns
  <key>_property.geojson  the boundary

Runs on GitHub Actions (which can reach the data). No invented numbers.
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

# ---- properties ------------------------------------------------------------
PROPERTIES = {
    "plevna":    {"name": "Plevna Downs",     "lots": ["6SP196201"]},
    "tilbooroo": {"name": "Tilbooroo Station", "lots": ["1336PH1989"]},
}
CELL_M = int(os.environ.get("CELL_M", "1000"))
OUT    = "app_data"; os.makedirs(OUT, exist_ok=True)
ALBERS = 3577
END    = date.today()
# four 7-day windows ending today: WEEKS[0] = (END-7, END) = this week
WEEKS  = [(END - timedelta(days=7*(i+1)), END - timedelta(days=7*i)) for i in range(4)]
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
    poly = prop_al.union_all() if hasattr(prop_al, "union_all") else prop_al.unary_union
    g = g[g.intersects(poly)].copy()
    g["geometry"] = g.intersection(poly)
    g = g[~g.is_empty].reset_index(drop=True)
    g["cell_id"] = [f"C{i:03d}" for i in range(len(g))]
    return g

def zmean(grid, href, crs, nodata):
    gg = grid.to_crs(crs)
    st = zonal_stats(gg, f"/vsicurl/{to_https(href)}", stats=["mean"], nodata=nodata, all_touched=True)
    return [None if s["mean"] is None else s["mean"] for s in st]

dea = Client.open("https://explorer.dea.ga.gov.au/stac")
mpc = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)

for key, meta in PROPERTIES.items():
    try:
        log("===", key, meta["lots"])
        prop = gpd.GeoDataFrame(geometry=[cadastre(meta["lots"])], crs=4326)
        prop.to_file(f"{OUT}/{key}_property.geojson", driver="GeoJSON")
        prop_al = prop.to_crs(ALBERS)
        grid = make_grid(prop_al)
        bbox = [float(v) for v in prop.total_bounds]
        log(f"{key}: {len(grid)} cells, {prop_al.area.iloc[0]/1e4:,.0f} ha")

        for i, (ws, we) in enumerate(WEEKS):
            tag = f"w{i}"
            # --- cover: latest DEA FC scene in the window ---
            try:
                items = list(dea.search(collections=["ga_ls_fc_3"], bbox=bbox,
                                        datetime=f"{ws}/{we}").items())
                if items:
                    it = max(items, key=lambda x: x.datetime)
                    crs = rasterio.open(f"/vsicurl/{to_https(it.assets['bs'].href)}").crs
                    bs = zmean(grid, it.assets["bs"].href, crs, 255)
                    pv = zmean(grid, it.assets["pv"].href, crs, 255)
                    grid[f"cover_{tag}"]     = [None if b is None else round(100 - b) for b in bs]
                    grid[f"green_{tag}"]     = [None if p is None else round(p) for p in pv]
                    grid[f"coverdate_{tag}"] = str(it.datetime.date())
                else:
                    grid[f"cover_{tag}"] = None; grid[f"green_{tag}"] = None; grid[f"coverdate_{tag}"] = None
            except Exception as e:
                log(f"{key} {tag} cover failed:", e)
                grid[f"cover_{tag}"] = None; grid[f"green_{tag}"] = None; grid[f"coverdate_{tag}"] = None
            # --- temp: latest Landsat scene in the window ---
            try:
                items = list(mpc.search(collections=["landsat-c2-l2"], bbox=bbox,
                             datetime=f"{ws}/{we}", query={"eo:cloud_cover": {"lt": 60}}).items())
                if items:
                    it = max(items, key=lambda x: x.datetime)
                    crs = rasterio.open(f"/vsicurl/{it.assets['lwir11'].href}").crs
                    dn = zmean(grid, it.assets["lwir11"].href, crs, 0)
                    grid[f"temp_{tag}"]     = [None if d is None else round(d*0.00341802 + 149.0 - 273.15, 1) for d in dn]
                    grid[f"tempdate_{tag}"] = str(it.datetime.date())
                else:
                    grid[f"temp_{tag}"] = None; grid[f"tempdate_{tag}"] = None
            except Exception as e:
                log(f"{key} {tag} temp failed:", e)
                grid[f"temp_{tag}"] = None; grid[f"tempdate_{tag}"] = None
            # --- rain: SILO daily totalled over the window (EPSG:4326) ---
            try:
                total = None; days = 0; d = ws
                while d < we:
                    u = f"{SILO_DAILY}{d.year}/{d.year}{d.month:02d}{d.day:02d}.daily_rain.tif"
                    try:
                        vals = zmean(grid, u, WGS84, None)
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
                log(f"{key} {tag} rain failed:", e)
                grid[f"rain_{tag}"] = None

        grid.to_crs(4326).to_file(f"{OUT}/{key}_grid.geojson", driver="GeoJSON")
        log(f"{key}: wrote grid, cols:", [c for c in grid.columns if c != "geometry"])
    except Exception as e:
        log(key, "FAILED:", e); traceback.print_exc()

# week labels for the viewer
import json
labels = {f"w{i}": (("This week" if i == 0 else f"{i} week{'s' if i>1 else ''} ago")
                    + f"  ({ws:%d %b}-{we:%d %b})") for i, (ws, we) in enumerate(WEEKS)}
json.dump({"weeks": labels, "built": str(END)}, open(f"{OUT}/meta.json", "w"))
log("Done.", labels)
