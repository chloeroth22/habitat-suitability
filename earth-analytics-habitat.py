# Import packages
import os
import pathlib
import zipfile
from urllib.parse import urljoin

import earthaccess  # Access NASA data from the cloud
import earthpy as et
import geopandas as gpd  # Geospatial data handling
import matplotlib.pyplot as plt  # Overlay raster and vector data
import numpy as np  # NumPy for array manipulation
import pandas as pd  # Group and aggregate data
import requests  # HTTP requests handling
import rioxarray as rxr  # Work with raster data
from rioxarray.merge import merge_arrays  # Merge rasters
from bs4 import BeautifulSoup  # BeautifulSoup for parsing HTML
import netCDF4  # Work with netCDF files
import xarray as xr  # Work with multi-dimensional arrays
import cartopy.crs as ccrs  # Cartopy for map projections
import seaborn as sns  # Statistical data visualization
import rasterio as rio  # Raster I/O
import regionmask  # Masking regions in geospatial data
from rasterio.crs import CRS  # Raster CRS handling
import xrspatial as xrs

# Create data directory in the home folder
data_dir = os.path.join(
        pathlib.Path.home(),
        'earth-analytics',
        'data',
        'habitat',
)
os.makedirs(data_dir, exist_ok=True)

# Set up grassland  boundary URL
grassland_url = (
        "https://data.fs.usda.gov/geodata/edw/edw_resources/shp/"
        "S_USA.NationalGrassland.zip")

# Set up path to save the grassland data on your machine
grassland_dir = os.path.join(data_dir, 'grassland')
os.makedirs(grassland_dir, exist_ok=True)
grassland_path = os.path.join(grassland_dir, 'grassland.shp')
grassland_north_path = os.path.join(grassland_dir, 'grassland_north.shp')
grassland_south_path = os.path.join(grassland_dir, 'grassland_south.shp')

# Only download grassland data once
if not os.path.exists(grassland_path):
    grassland_gdf = gpd.read_file(grassland_url)
    grassland_gdf.to_file(grassland_path)
    print(f"Downloading {grassland_url}...")

else:
    print(f"File {grassland_path} already exists. Skipping download.")

# Open up grasslands ande create shp for North, Sheyenne,
# and South, Lyndon B. Johnson, National Grasslands
grassland_gdf = gpd.read_file(grassland_path)
grassland_gdf.iloc[[3]].to_file(grassland_north_path)
grassland_north_gdf = gpd.read_file(grassland_north_path)
grassland_gdf.iloc[[10]].to_file(grassland_south_path)
grassland_south_gdf = gpd.read_file(grassland_south_path)

# Reproject grassland gdfs into WGS84 and UTM Zone 14
grassland_north_gdf_UTM = grassland_north_gdf.to_crs('EPSG:32614')
grassland_south_gdf_UTM = grassland_south_gdf.to_crs('EPSG:32614')
grassland_north_gdf_WGS = grassland_north_gdf.to_crs('EPSG:4326')
grassland_south_gdf_WGS = grassland_south_gdf.to_crs('EPSG:4326')

# Set up POLARIS url
polaris_url = ("http://hydrology.cee.duke.edu/POLARIS/PROPERTIES/"
               "v1.0/ph/mean/5_15/")

# Set up path to save POLARIS data to machine
polaris_dir = os.path.join(data_dir, 'polaris')
os.makedirs(polaris_dir, exist_ok=True)

# Send an HHTP GET request for HTML content of the page
response = requests.get(polaris_url)
response.raise_for_status()  # Raise am exeption if the request failed

# Parse the HTML
polaris_parse = BeautifulSoup(response.text, 'html.parser')

# Find all the ancho tags with href attributes
polaris_links = polaris_parse.find_all('a', href=True)

# Define latitude and longitude search patterns
lat_keywords = ['lat4647', 'lat3334']
lon_keywords = ['lon-98-97', 'lon-97-96']

# Filter out the links that end with '.tif' and contain both a latitude
# and a longitude pattern
polaris_tif_links = [
    urljoin(polaris_url, link['href'])
    for link in polaris_links
    if link['href'].endswith('.tif')  # Only .tif files
    and any(lat in link['href'] for lat in lat_keywords)
    and any(lon in link['href'] for lon in lon_keywords)
]


# Funtion to download a .tif file
def download_tif(url, save_path):
    """
    Downloads a .tif file from the given URL and saves it to the specified
    path and skips the download if the file aready exists at the path.

    Parameters
    ----------
    url : str
        The URL from which the .tif file will be downloaded.
    save_path : str
        The local path where the downloaded .tif file will be saved.

    Raises
    ------
    requests.exceptions.RequestException
        If an error occurs during the HTTP request a `RequestException` will be
        raised.

    """
    if os.path.exists(save_path):
        print(f"File {save_path} already exists. Skipping download.")
        return
    print(f"Downloading {url}...")
    response = requests.get(url)
    response.raise_for_status()

    # Write the content of the file to the specified path
    with open(save_path, 'wb') as f:
        f.write(response.content)


# Initalize POLARIS lists
polaris_files_list_north = []
polaris_files_list_south = []

# Download all .tif files from POLARIS website
for tif_url in polaris_tif_links:
    polaris_tif_filename = os.path.basename(tif_url)

    # Create the full file path where the file will be saved
    polaris_save_path = os.path.join(polaris_dir, polaris_tif_filename)

    # Create list of all files in POLARIS folder, sorting by North and South
    if (lat_keywords[1] in polaris_save_path and lon_keywords[0]
        in polaris_save_path):
        polaris_files_list_south.append(polaris_save_path)
    if lat_keywords[0] in polaris_save_path:
        polaris_files_list_north.append(polaris_save_path)

    # Download and save the .tif file
    download_tif(tif_url, polaris_save_path)

# Read rasters into xarray
polaris_north_rasters = [
        rxr.open_rasterio(file, masked=True).squeeze()
        for file in polaris_files_list_north
]
polaris_south_rasters = [
        rxr.open_rasterio(file, masked=True).squeeze()
        for file in polaris_files_list_south
]

# Crop each raster to the grasslands .shp
polaris_north_cropped_rasters = [
        raster.rio.clip(
            grassland_north_gdf.geometry.values,
            grassland_north_gdf.crs)
        for raster in polaris_north_rasters
]
polaris_south_cropped_rasters = [
        raster.rio.clip(
            grassland_south_gdf.geometry.values,
            grassland_south_gdf.crs)
        for raster in polaris_south_rasters
]

# Merge the cropped rasters
polaris_north_merged_raster = merge_arrays(polaris_north_cropped_rasters)
polaris_south_merged_raster = merge_arrays(polaris_south_cropped_rasters)

# Plot regions
fig, ax = plt.subplots(figsize=(10, 8))
pc = polaris_north_merged_raster.plot(ax=ax, cmap='viridis')
ax.set_title('Soil pH — Sheyenne National Grassland')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.show()

fig, ax = plt.subplots(figsize=(10, 8))
polaris_south_merged_raster.plot(ax=ax, cmap='viridis')
ax.set_title('Soil pH — Lyndon B. Johnson National Grassland')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.show()

# Set up earthaccess directories
SRTM_dir = os.path.join(data_dir, 'SRTM')
os.makedirs(SRTM_dir, exist_ok=True)
SRTM_north_dir = os.path.join(SRTM_dir, 'SRTM_north')
os.makedirs(SRTM_north_dir, exist_ok=True)
SRTM_south_dir = os.path.join(SRTM_dir, 'SRTM_south')
os.makedirs(SRTM_south_dir, exist_ok=True)

# Set up earthaccess conncection
earthaccess.login(strategy="interactive", persist=True)

# Initalize lists to store paths of unzipped files
SRTM_north_files_list = []
SRTM_south_files_list = []


# Function to check if data has been downloaded
def is_data_downloaded(directory):
    """
    Checks if any files are present in the specified directory, indicating
    that data has been downloaded.

    Parameters
    ----------
    directory : str
        The path to the directory that is being checked for downloaded files.

    Returns
    -------
    bool
        Returns `True` if the directory contains any files, indicating that
        data has been downloaded. Returns `False` if the directory is empty.
    """
    return len(os.listdir(directory)) > 0


# Search and download earthaccess data for North boundary
if not is_data_downloaded(SRTM_north_dir):
    SRTM_north_results = earthaccess.search_data(
        doi="10.5067/MEaSUREs/SRTM/SRTMGL1.003",
        bounding_box=tuple(grassland_north_gdf_WGS.total_bounds),
        )
    SRTM_north_files = earthaccess.download(SRTM_north_results, SRTM_north_dir)
    # Unzip then delete .zip files
    for file in SRTM_north_files:
        with zipfile.ZipFile(file, 'r') as zip_ref:
            zip_ref.extractall(SRTM_north_dir)

        # After extracting, add the unzipped file(s) to the list
        for extracted_file in zip_ref.namelist():
            extracted_file_path = os.path.join(SRTM_north_dir, extracted_file)
            SRTM_north_files_list.append(extracted_file_path)

        # Remove the zip file after extraction
        if file.endswith('.zip'):
            os.remove(file)
else:
    # If the data is already downloaded, populate the file list directly
    SRTM_north_files_list = [
            os.path.join(SRTM_north_dir, f) for f in os.listdir(SRTM_north_dir)
            if f.endswith('.hgt')
    ]
    print(f"Files {SRTM_north_files_list} already exists. Skipping download.")

# Search and download earthaccess data for South boundary
if not is_data_downloaded(SRTM_south_dir):
    SRTM_south_results = earthaccess.search_data(
        doi="10.5067/MEaSUREs/SRTM/SRTMGL1.003",
        bounding_box=tuple(grassland_south_gdf_WGS.total_bounds),
    )
    SRTM_south_files = earthaccess.download(SRTM_south_results, SRTM_south_dir)

    # Unzip then delete .zip files
    for file in SRTM_south_files:
        with zipfile.ZipFile(file, 'r') as zip_ref:
            zip_ref.extractall(SRTM_south_dir)

        # After extracting, add the unzipped file(s) to the list
        for extracted_file in zip_ref.namelist():
            extracted_file_path = os.path.join(SRTM_south_dir, extracted_file)
            SRTM_south_files_list.append(extracted_file_path)

        # Remove the zip file after extraction
        if file.endswith('.zip'):
            os.remove(file)
else:
    # If the data is already downloaded, populate the file list directly
    SRTM_south_files_list = [
            os.path.join(SRTM_south_dir, f) for f in os.listdir(SRTM_south_dir)
            if f.endswith('.hgt')
    ]
    print(f"Files {SRTM_south_files_list} already exists. Skipping download.")

# Read rasters into xarray
SRTM_south_rasters = [
        rxr.open_rasterio(file, masked=True).squeeze()
        for file in SRTM_south_files_list
]
SRTM_north_rasters = [
        rxr.open_rasterio(file, masked=True).squeeze()
        for file in SRTM_north_files_list
        ]

# Crop each raster to the grasslands .shp
SRTM_south_cropped_rasters = [
    raster.rio.clip(grassland_south_gdf.geometry.values,
                    grassland_south_gdf.crs)
    for raster in SRTM_south_rasters
]
SRTM_north_cropped_rasters = [
    raster.rio.clip(grassland_north_gdf.geometry.values,
                    grassland_north_gdf.crs)
    for raster in SRTM_north_rasters
]

# Merge the cropped rasters
SRTM_south_merged_raster = merge_arrays(SRTM_south_cropped_rasters)
SRTM_north_merged_raster = merge_arrays(SRTM_north_cropped_rasters)

# Plot regions
fig, ax = plt.subplots(figsize=(10, 8))
pc = SRTM_north_merged_raster.plot(ax=ax, cmap='viridis')
ax.set_title('Elevation (meters) — Sheyenne National Grassland')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.show()

fig, ax = plt.subplots(figsize=(10, 8))
SRTM_south_merged_raster.plot(ax=ax, cmap='viridis')
ax.set_title('Elevation (meters)— Lyndon B. Johnson National Grassland')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.show()

# Reproject into UTM to calculate slope
crs_utm = CRS.from_string('EPSG:32614')
SRTM_north_merged_repro = SRTM_north_merged_raster.rio.reproject(crs_utm)
SRTM_south_merged_repro = SRTM_south_merged_raster.rio.reproject(crs_utm)

# Calculate slope
SRTM_south_slope = xrs.slope(SRTM_south_merged_repro)
SRTM_north_slope = xrs.slope(SRTM_north_merged_repro)

# Download MACA data
maca_url = ("http://thredds.northwestknowledge.net:8080/thredds/dodsC/"
            "agg_macav2metdata_pr_BNU-ESM_r1i1p1_historical_1950_2005"
            "_CONUS_monthly.nc")

# Open the data from the thredds server
try:
    precip = xr.open_dataset(maca_url)

except OSError as oe:
    print("Oops, it looks like the file that you are trying to connect to, "
          "{}, doesn't exist. Try to revisit your model options to ensure "
          "the data exist on the server.  ".format(maca_url))


# Function to calculate Lat and Lon of AOI
def get_aoi(shp, world=True):
    """
    Extracts the bounding box (latitude and longitude) of the area of interest
    (AOI) from a shapefile and returns it as a dictionary.

    Parameters
    ----------
    shp : geopandas.GeoDataFrame
        A GeoDataFrame containing the geometry of the AOI.
        The `total_bounds` attribute of this GeoDataFrame is used to determine
        the lat and lon bounds.

    world : bool, optional, default=True
        If `True`, the function will adjust the longitude values by adding 360
        to both the min and max longitudes to account for global coordinates.
        If `False`, no adjustment is made.

    Returns
    -------
    dict
        A dictionary containing the latitude and longitude bounds of the AOI:
        - "lat" : list of float
            A list of two floats representing the min and max latitudes.
        - "lon" : list of float
            A list of two floats representing the min and max longitudes.
    """
    lon_lat = {}
    # Get lat min, max
    aoi_lat = [float(shp.total_bounds[1]), float(shp.total_bounds[3])]
    aoi_lon = [float(shp.total_bounds[0]), float(shp.total_bounds[2])]

    if world:
        aoi_lon[0] = aoi_lon[0] + 360
        aoi_lon[1] = aoi_lon[1] + 360
    lon_lat["lon"] = aoi_lon
    lon_lat["lat"] = aoi_lat
    return lon_lat


# Calculate Southern boundary
south_bounds = get_aoi(grassland_south_gdf_WGS)

# Calculate Northern boundary
north_bounds = get_aoi(grassland_north_gdf_WGS)


# Funtion to slice, mask, and average the precipitation data
def conus(dataset, start_date, end_date, bounds, gdf):
    """
    Slices, masks, and averages precipitation data over a specified time period
    and spatial region, then calculates the average precipitation.

    Parameters
    ----------
    dataset : xarray.Dataset
        Precipitation dataset with time, lat, and long dimensions.
        Assumed that the dataset contains a variable called "precipitation".

    start_date : str or datetime-like
        The start date for slicing the data, formatted as 'YYYY-MM-DD'

    end_date : str or datetime-like
        The end date for slicing the data, formatted as 'YYYY-MM-DD'

    bounds : dict
        A dictionary containing the spatial bounds for slicing the data:
        - "lat" : list of float
            A list of two floats specifying the min and max latitude bounds.
        - "lon" : list of float
            A list of two floats specifying the min and max longitude bounds.

    gdf : geopandas.GeoDataFrame
        A GeoDataFrame containing geometries used to mask precipitation data.

    Returns
    -------
    xarray.DataArray
        An xarray DataArray containing the average precipitation over the time
        period and spatial extent, masked, and reprojected to EPSG:4326 CRS.
    """
    # Slice the data by time and spatial extent
    conus_data = dataset["precipitation"].sel(
        time=slice(start_date, end_date),
        lon=slice(bounds["lon"][0], bounds["lon"][1]),
        lat=slice(bounds["lat"][0], bounds["lat"][1]))

    # Create mask for region
    mask = regionmask.mask_3D_geopandas(gdf,
                                        conus_data.lon,
                                        conus_data.lat)
    # Apply the mask to the data
    conus_data_masked = conus_data.where(mask)

    # Calculate average precipitation over time
    conus_data_summary = conus_data_masked.mean(dim="time", skipna=True)
    conus_data_summary = conus_data_summary.transpose('region', 'lat', 'lon')
    conus_data_summary.rio.write_crs("EPSG:4326", inplace=True)
    return conus_data_summary


# Run conus funtion on each region and time frame
south_1950 = conus(precip, "1950-01-15", "1954-12-15",
                   south_bounds, grassland_south_gdf_WGS)
south_2000 = conus(precip, "2000-01-15", "2004-12-15", south_bounds,
                   grassland_south_gdf_WGS)
north_1950 = conus(precip, "1950-01-15", "1954-12-15", north_bounds,
                   grassland_north_gdf_WGS)
north_2000 = conus(precip, "2000-01-15", "2004-12-15", north_bounds,
                   grassland_north_gdf_WGS)
# Plot regions
fig, ax = plt.subplots(figsize=(10, 8))
pc = north_1950.plot(ax=ax, cmap='viridis')
ax.set_title('Average Precipitation (mm) 1950 to 1954 — '
             'Sheyenne National Grassland')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.show()

fig, ax = plt.subplots(figsize=(10, 8))
pc = north_2000.plot(ax=ax, cmap='viridis')
ax.set_title('Average Precipitation (mm) 2000 to 2004 — '
             'Sheyenne National Grassland')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.show()

fig, ax = plt.subplots(figsize=(10, 8))
pc = south_1950.plot(ax=ax, cmap='viridis')
ax.set_title('Average Precipitation (mm) 1950 to 1954 — '
             'Lyndon B. Johnson National Grassland')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.show()

fig, ax = plt.subplots(figsize=(10, 8))
pc = south_2000.plot(ax=ax, cmap='viridis')
ax.set_title('Average Precipitation (mm) 2000 to 2004 — '
             'Lyndon B. Johnson National Grassland')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
plt.show()


# Function to harmonize data
def match(xds, xds_match):
    """
    Harmonizes the coordinate system and dimensions of two xarray datasets.

    Parameters
    ----------
    xds : xarray.DataArray or xarray.Dataset
        The first xarray dataset or data array to be harmonized. It is the reference dataset 
        that will be aligned with the second dataset.

    xds_match : xarray.DataArray or xarray.Dataset
        The second xarray dataset or data array that will be reprojected to match the coordinate 
        system of `xds`. It is assumed to have compatible spatial dimensions (`x` and `y`).

    Returns
    -------
    xarray.DataArray or xarray.Dataset
        A new xarray object with the same data as `xds`, but with coordinates aligned to match 
        the `xds_match` dataset's coordinates (using re-projection). The coordinates for `x` 
        and `y` are reassigned based on `xds_match`.
    """
    xds_repr_match = xds.rio.reproject_match(xds_match)
    xds_repr_match = xds_repr_match.assign_coords({
        "x": xds_match.x,
        "y": xds_match.y,
    })
    return xds_repr_match


# Harmonize data
polaris_n = match(polaris_north_merged_raster, SRTM_north_slope)
n_1950 = match(north_1950, SRTM_north_slope)
n_2000 = match(north_2000, SRTM_north_slope)
polaris_s = match(polaris_south_merged_raster, SRTM_south_slope)
s_1950 = match(south_1950, SRTM_south_slope)
s_2000 = match(south_2000, SRTM_south_slope)


# Create a fuzzy logic habitat suitability model

def linear(x, slope, intercept):
    """
    Linear funtion to calculate value between 0 and 1 for values
    outside the optimal range.
    Args:
    -x: The value to apply the linear funtion to.
    -slope: The slope of the line
    -intercept: The intercept of the line.

    Returns:
    -Scaled value based on the linear function.
    """
    return slope * x + intercept


def transform_raster(raster, min_dealbreaker, min_optimal, max_optimal,
                     max_dealbreaker):
    """
    Applies a linear transformation to a raster dataset based on specified threshold values.

    The function normalizes the values in the raster by applying a linear transformation
    to values that fall within specified ranges. The transformation is defined using
    two linear functions: one for values below the optimal range (min) and one for values
    above the optimal range (max). Values outside of the defined ranges are assigned boundary 
    values (0 or 1).

    Parameters
    ----------
    raster : xarray.DataArray
        The input raster dataset to be transformed.

    min_dealbreaker : float
        The minimum value threshold below which all raster values are mapped to 0.

    min_optimal : float
        The minimum value threshold above which raster values will be linearly transformed 
        using the `slope_min` and `intercept_min` values.

    max_optimal : float
        The maximum value threshold below which raster values will be linearly transformed 
        using the `slope_max` and `intercept_max` values.

    max_dealbreaker : float
        The maximum value threshold above which all raster values are mapped to 0.

    Returns
    -------
    xarray.DataArray
        A new xarray DataArray with the transformed raster values, where values outside the 
        specified ranges are clipped to 0 or 1, and values within the optimal range are 
        linearly transformed.

    Notes
    -----
    - Values below `min_dealbreaker` are mapped to 0, and values above `max_dealbreaker`
      are also mapped to 0.
    - Values between `min_dealbreaker` and `min_optimal` are linearly scaled from 0 to 1.
    - Values between `max_optimal` and `max_dealbreaker` are linearly scaled from 1 to 0.
    - Values between `min_optimal` and `max_optimal` are mapped to 1.
    """
    if not min_optimal == 0:
        slope_min = (1-0)/(min_optimal-min_dealbreaker)
        intercept_min = 0-(slope_min*min_dealbreaker)
    slope_max = (0-1)/(max_dealbreaker-max_optimal)
    intercept_max = 1-(slope_max*max_optimal)

    def transform(value):
        if np.isnan(value):
            return np.nan
        if value < min_dealbreaker:
            return 0
        if value < min_optimal:
            return linear(value, slope_min, intercept_min)
        if value > max_dealbreaker:
            return 0
        if value > max_optimal:
            return linear(value, slope_max, intercept_max)
        else:
            return 1
    # Apply the normalization function to each raster value
    transformed_raster = raster.copy()
    transformed_raster.values = np.vectorize(transform)(raster.values)

    return transformed_raster


# S. Nutans does best at a pH of pH range of 4.8 to 8.0
# ideally between 5 and 6.
polaris_n_tran = transform_raster(polaris_n, 4.8, 5, 6, 8)
polaris_s_tran = transform_raster(polaris_s, 4.8, 5, 6, 8)

# S. Nutans does best at a monthly rainfall of 42mm to 64mm
# (or annual of 28 cm to 114 cm)
n_1950_tran = transform_raster(n_1950, 23, 42, 64, 95)
n_2000_tran = transform_raster(n_2000, 23, 42, 64, 95)
s_1950_tran = transform_raster(s_1950, 23, 42, 64, 95)
s_2000_tran = transform_raster(s_2000, 23, 42, 64, 95)

# S. Nutans does best at low angle slopes
SRTM_s_tran = transform_raster(SRTM_south_slope, 0, 0, 30, 90)
SRTM_n_tran = transform_raster(SRTM_north_slope, 0, 0, 30, 90)

# Multipy raters together
north_optimized_1950 = polaris_n_tran * n_1950_tran * SRTM_n_tran
north_optimized_2000 = polaris_n_tran * n_2000_tran * SRTM_n_tran
south_optimized_1950 = polaris_s_tran * s_1950_tran * SRTM_s_tran
south_optimized_2000 = polaris_s_tran * s_2000_tran * SRTM_s_tran

fig, ax = plt.subplots(figsize=(10, 8))
pc = north_optimized_1950.plot(ax=ax, cmap='viridis')
ax.set_title('Habitat Suitability 1950 to 1954 — '
             'Sheyenne National Grassland')
ax.set_xlabel('Meters')
ax.set_ylabel('Meters')
plt.show()

fig, ax = plt.subplots(figsize=(10, 8))
pc = north_optimized_2000.plot(ax=ax, cmap='viridis')
ax.set_title('Habitat Suitability 2000 to 2004 — '
             'Sheyenne National Grassland')
ax.set_xlabel('Meters')
ax.set_ylabel('Meters')
plt.show()

fig, ax = plt.subplots(figsize=(10, 8))
pc = south_optimized_1950.plot(ax=ax, cmap='viridis')
ax.set_title('Habitat Suitability 1950 to 1954 — '
             'Lyndon B. Johnson National Grassland')
ax.set_xlabel('Meters')
ax.set_ylabel('Meters')
plt.show()

fig, ax = plt.subplots(figsize=(10, 8))
pc = south_optimized_2000.plot(ax=ax, cmap='viridis')
ax.set_title('Habitat Suitability 2000 to 2004 — '
             'Lyndon B. Johnson National Grassland')
ax.set_xlabel('Meters')
ax.set_ylabel('Meters')
plt.show()
