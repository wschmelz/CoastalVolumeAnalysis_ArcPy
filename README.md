# ArcPy Tool for Coastal Geomorphological Change Analysis

## Overview
The 3D_Analyst.py script programatically processes multiple temporally distinct topographical survey datasets and calculates the volumetric change occurring between surveys throughout a survey site using ArcGIS and ArcPy. The script processes survey data stored within a shapefile to generate Digital Elevation Models (DEMs) within a defined three-dimensional spatial domain, and makes the volumetric change calculations. The outputs of the script are data tables that record spatiotemporal calculations of volumetric change, raster datasets that represent the spatial distribution of erosion and deposition between surveys, and images that visualize the spatiotemporal changes. 

The Python based analytical tool we provide performs the data analysis tasks defined by the “Coastal Landform Elevation Models” protocol of the Northeast Coastal and Barrier Network of the National Park Service. Please cite this protocol in the documentation of any application of this software: 

Psuty, N. P., W. J. Schmelz, and A. Habeck. 2018. Northeast Coastal and Barrier Network 
geomorphological monitoring protocol: Part III – coastal landform elevation models. Natural
Resource Report NPS/NCBN/NRR—2018/1712. National Park Service, Fort Collins, Colorado.

**GitHub Repository**: [CoastalVolumeAnalysis_ArcPy](https://github.com/wschmelz/CoastalVolumeAnalysis_ArcPy.git)

## Key methodological aspects
- DEM interpolation: Script generates DEMs from topographical survey data using Delaunay triangulation.
- Volumetric change calculation: Script generates difference rasters between surveys and uses these data to calculate volumetric change across the survey area.
- Calculation of intrasite spatial variability of voumetric change: Volumes of erosion and deposition are discretely calculated for each of any number of spatial compartments within the study area.

**Output Tables and Maps**:
- Tables containing metrics of volumetric changes.
- Raster and feature class spatial data that portray spatial distribution of erosion and deposition through time.
- ArcGIS map projects with images portraying the volumetric change that were generated using the survey data.

## Repository Structure
```
/01_Benchmarks/           # Directory to store benchmark control points for surveys (not needed to run script)
/02_RawData/              # Directory to store raw input survey data (not needed to run script)
/03_Shapefiles/           # Input shapefiles that contain topographical survey data
/04_TINs/                 # Directory that stores script generated Triangular Irregular Networks (TINs)
/05_ArcMap/               # Directory that stores ArcGIS project files for visualization, including the template ArcGIS project that the script updates based on data
/06_Rasters/              # Directory that stores raster datasets genetrated from the TINs
/07_Differences/          # Directory that stores difference rasters (elevation changes)
/08_Tables/               # Tables that store volumetric change calculations
/09_Compartments/         # Directory that stores shapefile with compartments for use in volumetric change analysis
/10_Vectors/              # Directory that stores geodatabase with feature classes that portray spatial variability in volumetric changes alongshore within the site
/11_Layers/               # Directory that stores layer files for ArcGIS template and visualization
/12_TotalDifferences/     # Endpoint difference rasters
/13_TotalTables/          # Endpoint volumetric change  tables
/14_TotalVectors/         # Endpoint vector feature classes
/15_Images/               # Map images portraying volumetric change
3D_Analyst.py             # Main script
3D_Analyst_params.txt     # Analysis and map visualization parameter file
Compartment_params.txt    # Compartment/baseline parameter file
```

## Required software
- ArcGIS Pro (with Spatial Analyst and 3D Analyst extensions)
- Python 3.x with ArcPy installed

## Data
Input survey data must be provided as shapefiles named in the format: `[Code]_3D_[YYYYMMDD].shp` (e.g., `CZo_3D_20140131.shp`).

## Configuration Files
- `3D_Analyst_params.txt`: Specifies analysis parameters.
- `Compartment_params.txt`: Defines spatial compartment geometry.

## Setup and Usage

### Step 1: Prepare Data
- Create shapefiles from survey data with parameters defined in Psuty et al. (2018) SOP6, and place these survey shapefiles into the `/03_Shapefiles/` directory.
- Create data items like the base surface and compartment shapefiles that are described in Psuty et al. (2018) SOP7 
- Ensure shapefiles adhere to the naming format and include the required attributes `[ID,X,Y,Z,CODE]`.

### Step 2: Configure Parameters
Edit the following files based on the study area and analysis requirements:

**3D_Analyst_params.txt: Example**
```
vectormultiplier = .08          # multiplier for visualization vectors
compbreak1 = 15                 # spatial breakpoint 1
compbreak2 = 30                 # spatial breakpoint 2
startdate = "20140131"          # initial survey date
enddate = "20160121"            # final survey date
marker_labels = 250,500,1000    # volume in cubic meters that a image grid label will be created for
cellsize = 1.0                  # raster cell size
```

**Compartment_params.txt: Example**
```
easting = 586681        # x-coordinate of baseline origin vertex
northing = 4474592      # y-coordinate of baseline origin vertex
baselineangle = 97.5    # angle of baseline (cw degrees relative to east)
compdirection = "left"  # direction of shoreline from the basline
distance = 30           # alongshore distance for each compartment
length = 700            # max cross-shore length of compartments
compartments = 45       # number of compartments in the study area
elev_threshold = -0.335 # elevation threshold (all topography below this threshold will not be considered in volumetric change analysis)
```

### Step 3: Run the Script

```
python 3D_Analyst.py
```

## Output
The script will produce:

- DEM and difference rasters located in `/06_Rasters/`, `/07_Differences/`, and `/12_TotalDifferences/`.
- Summary tables in excel format in `/08_Tables/` and `/13_TotalTables/`.
- Feature classes with vectors that reflect volumetric change at points alongshore in `/10_Vectors/` and `/14_TotalVectors/`.
- Map images of spatial distribution of volumetric change between each survey and across the entire study period in `/15_Images/`.
