#!/usr/bin/env python3
"""
build_grid.py  -  real gridded paddock condition for a property.

Runs on GitHub Actions (which can reach the data). No invented numbers:
every cell value is the mean of the actual satellite pixels inside it.

Writes into app_data/:
  property.geojson    the property boundary (from the QLD cadastre)
  cover_grid.geojson  grid cells with real ground-cover %, green %, temperature
  rain_7day.txt       one 7-day rainfall figure (SILO ~5 km), if SILO_EMAIL set

Config comes from environment variables (set in the workflow):
  LOTS        comma-separated lot/plans, e.g. "6SP196201,4WLA10"
  CELL_M      grid cell size in metres (default 1000)
  SILO_EMAIL  your email = SILO's free API key (optional)
"""
import os, json, traceback
from datetime import date, timedelta

import requests
import geopandas as gpd
from shapely.geometry import shape, box
from shapely.ops import unary_union
from pystac_client import Client
import planetary_computer as pc
import rasterio                      # noqa: F401 (ensures GDAL/curl present)
from rasterstats import zonal_stats

LOTS   = [s.strip().upper() for s in os.environ.get("LOTS", "6SP196201").split(",") if s.strip()]
CELL_M = int(os.environ.get("CELL_M", "1000"))
EMAIL  = os.environ.get("SILO_EMAIL", "").strip()
OUT    = "app_data"; os.makedirs(OUT, exist_ok=True)
ALBERS = 3577                        # Australian Albers: equal-area, DEA's own CRS
END    = date.today(); START = END - timedelta(days=120)

def log(*a): print(*a, flush=True)

# -------------------------------------------------- 1. property from cadastre
log(f"Lots: {LOTS}")
where = "lotplan IN (%s)" % ",".join(f"'{l}'" for l in LOTS)
cad = "https://spatial-gis.information.qld.gov.au/arcgis/rest/services/Basemaps/FoundationData/FeatureServer/2/query"
r = requests.get(cad, params=dict(where=where, outFields="lotplan",
                 returnGeometry="true", outSR=4326, f="geojson"), timeout=90)
r.raise_for_status()
feats = r.json().get("features", [])
if not feats:
    raise SystemExit(f"No parcels matched {LOTS}. Check the lot/plan(s).")
prop = unary_union([shape(f["geometry"]) for f in feats])
prop_gdf = gpd.GeoDataFrame(geometry=[prop], crs=4326)
prop_gdf.to_file(f"{OUT}/property.geojson", driver="GeoJSON")
bbox = [float(x) for x in prop_gdf.total_bounds]        # 4326 minx,miny,maxx,maxy
log(f"Property: {prop_gdf.to_crs(ALBERS).area.iloc[0]/1e4:,.0f} ha from {len(feats)} part(s)")

# -------------------------------------------------- 2. grid (Albers metres)
prop_al = prop_gdf.to_crs(ALBERS)
minx, miny, maxx, maxy = prop_al.total_bounds
cells = []
x = minx
while x < maxx:
    y = miny
    while y < maxy:
        cells.append(box(x, y, x + CELL_M, y + CELL_M)); y += CELL_M
    x += CELL_M
grid = gpd.GeoDataFrame(geometry=cells, crs=ALBERS)
poly = prop_al.union_all() if hasattr(prop_al, "union_all") else prop_al.unary_union
grid = grid[grid.intersects(poly)].copy()
grid["geometry"] = grid.intersection(poly)              # trim cells to the boundary
grid = grid[~grid.is_empty].reset_index(drop=True)
grid["cell_id"] = [f"C{i:03d}" for i in range(len(grid))]
log(f"Grid: {len(grid)} cells of ~{CELL_M*CELL_M/1e4:.0f} ha")

def zonal_mean(href, crs, nodata):
    """Mean raster value per grid cell (grid reprojected to the raster's CRS)."""
    g = grid.to_crs(crs)
    stats = zonal_stats(g, f"/vsicurl/{href}", stats=["mean"], nodata=nodata, all_touched=True)
    return [None if s["mean"] is None else s["mean"] for s in stats]

# -------------------------------------------------- 3. DEA Fractional Cover -> cover %
try:
    dea = Client.open("https://explorer.dea.ga.gov.au/stac")
    items = list(dea.search(collections=["ga_ls_fc_3"], bbox=bbox,
                            datetime=f"{START}/{END}").items())
    it = max(items, key=lambda i: i.datetime)
    crs = rasterio.open(f"/vsicurl/{it.assets['bs'].href}").crs
    bs = zonal_mean(it.assets["bs"].href, crs, 255)     # bare soil %
    pv = zonal_mean(it.assets["pv"].href, crs, 255)     # green veg %
    grid["cover_pct"] = [None if b is None else round(100 - b) for b in bs]   # total ground cover
    grid["green_pct"] = [None if p is None else round(p) for p in pv]
    grid["fc_date"]   = str(it.datetime.date())
    log(f"Ground cover from DEA FC scene {it.datetime.date()}")
except Exception as e:
    log("FC step failed:", e); traceback.print_exc()

# -------------------------------------------------- 4. Landsat surface temperature (MPC)
try:
    mpc = Client.open("https://planetarycomputer.microsoft.com/api/stac/v1", modifier=pc.sign_inplace)
    items = list(mpc.search(collections=["landsat-c2-l2"], bbox=bbox,
                            datetime=f"{START}/{END}", query={"eo:cloud_cover": {"lt": 60}}).items())
    it = max(items, key=lambda i: i.datetime)
    crs = rasterio.open(f"/vsicurl/{it.assets['lwir11'].href}").crs
    dn = zonal_mean(it.assets["lwir11"].href, crs, 0)
    grid["temp_c"]    = [None if d is None else round(d*0.00341802 + 149.0 - 273.15, 1) for d in dn]
    grid["temp_date"] = str(it.datetime.date())
    log(f"Temperature from Landsat scene {it.datetime.date()}")
except Exception as e:
    log("Temperature step failed:", e); traceback.print_exc()

# -------------------------------------------------- 5. write grid (WGS84 for the map)
grid.to_crs(4326).to_file(f"{OUT}/cover_grid.geojson", driver="GeoJSON")
log("Wrote cover_grid.geojson with columns:", [c for c in grid.columns if c != "geometry"])

# -------------------------------------------------- 6. rainfall: one property figure (SILO ~5 km)
rain = "n/a"
if EMAIL:
    try:
        c = prop_gdf.to_crs(4326).geometry.iloc[0].centroid
        url = ("https://www.longpaddock.qld.gov.au/cgi-bin/silo/DataDrillDataset.php"
               f"?format=csv&start={(END-timedelta(days=7)):%Y%m%d}&finish={END:%Y%m%d}"
               f"&lat={c.y:.4f}&lon={c.x:.4f}&comment=R&username={EMAIL}&password=apirequest")
        rows = requests.get(url, timeout=60).text.strip().splitlines()
        hdr = rows[0].split(","); idx = hdr.index("daily_rain")
        rain = round(sum(float(r.split(",")[idx]) for r in rows[1:] if r.strip()), 1)
    except Exception as e:
        log("Rain step failed:", e)
open(f"{OUT}/rain_7day.txt", "w").write(f"Property 7-day rainfall (SILO ~5km): {rain} mm\n")
log("Done.")
