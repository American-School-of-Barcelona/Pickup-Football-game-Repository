import pandas as pd
import geopandas as gpd
import osmnx as ox
import numpy as np
import folium

csv_path = "football_pitches_strict.csv"
df = pd.read_csv(csv_path)
df = df.dropna(subset=["geo_epgs_4326_lat", "geo_epgs_4326_lon"])

gdf = gpd.GeoDataFrame(
    df,
    geometry=gpd.points_from_xy(df["geo_epgs_4326_lon"], df["geo_epgs_4326_lat"]),
    crs="EPSG:4326"
)

minx, miny, maxx, maxy = gdf.total_bounds
north, south, east, west = maxy, miny, maxx, minx

tags = {
    "leisure": "pitch",
    "sport": ["football", "soccer"]
}

osm = ox.geometries_from_bbox(north, south, east, west, tags)
osm_clean = osm[osm.geometry.notnull()].copy()
osm_clean["geom_point"] = osm_clean.geometry.centroid

osm_points = gpd.GeoDataFrame(
    osm_clean.drop(columns="geometry"),
    geometry="geom_point",
    crs=osm.crs
)

gdf_25831 = gdf.to_crs("EPSG:25831")
osm_25831 = osm_points.to_crs("EPSG:25831")

joined = gpd.sjoin_nearest(
    gdf_25831,
    osm_25831,
    how="left",
    distance_col="nearest_m"
)

DISTANCE_THRESHOLD = 40
joined["is_suspect"] = joined["nearest_m"] > DISTANCE_THRESHOLD

if "geo_epgs_25831_x" in joined.columns and "geo_epgs_25831_y" in joined.columns:
    dx = joined.geometry.x - joined["geo_epgs_25831_x"]
    dy = joined.geometry.y - joined["geo_epgs_25831_y"]
    joined["internal_offset_m"] = np.sqrt(dx**2 + dy**2)

joined_wgs84 = joined.to_crs("EPSG:4326")
osm_wgs84 = osm_25831.to_crs("EPSG:4326")

center_lat = joined_wgs84["geo_epgs_4326_lat"].mean()
center_lon = joined_wgs84["geo_epgs_4326_lon"].mean()

m = folium.Map(location=[center_lat, center_lon], zoom_start=12)

for _, row in joined_wgs84.iterrows():
    color = "red" if row["is_suspect"] else "green"
    tooltip = f"{row['name']} (nearest_m={row['nearest_m']:.1f})"
    folium.CircleMarker(
        location=[row.geometry.y, row.geometry.x],
        radius=5,
        color=color,
        fill=True,
        fill_opacity=0.9,
        tooltip=tooltip
    ).add_to(m)

for _, row in osm_wgs84.iterrows():
    geom = row.geometry
    if geom.geom_type in ["Polygon", "MultiPolygon"]:
        folium.GeoJson(geom, tooltip="OSM football pitch").add_to(m)
    else:
        folium.CircleMarker(
            location=[geom.y, geom.x],
            radius=3,
            color="blue",
            fill=True,
            fill_opacity=0.7,
            tooltip="OSM football pitch (point)"
        ).add_to(m)

joined_out = joined_wgs84.drop(columns="geometry")
joined_out.to_csv("football_pitches_osm_checked.csv", index=False)

output_html = "football_pitches_osm_check.html"
m.save(output_html)
