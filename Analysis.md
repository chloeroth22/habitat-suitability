## Habitat Suitability 
### Project Description  
For this project, I created a habitat suitability model for Sorghastrum nutans. S. nutans, known commonly as Indiangrass is native to North America, and is distributed from the east coast to Arizona, Utah, and the Rocky Mountains. 
It’s range has been moving northward over the last 50 years. 

Insert picture of grass and united states distribution

### Site Description
I choose two grasslands from the USFS National Grasslands to compare. To further investigate if S. nutans is moving northwards, I choose one National Grassland in the Northern United States and one in the Southern United States. 
The Northern grassland I picked was the Sheyenne National Grassland, located in North Dakota. Sheyenne National Grassland is 70,446 acres and located at 46 degrees North. The Southern Grassland that I picked was Lyndon B. Johnson 
National Grassland, located in Texas. Lyndon B. Johnson National Grassland is 20,309 acres and located at 33 degrees South.

### Data Sources
To determine habitat suitability I used data on soil pH, elevation (to calculate slope), and historical precipitation.
#### Soil pH
Soil pH was taken from POLARIS, which is a complete map of soil series probabilities that has been produced for the contiguous United States. POLARIS uses geospatial enviornmental data and a machine learning algorihm (DSMART-HPC) to
remap the Soil Survy Geographic (SSURGO) database. Data is available on the [POLARIS dataset website](http://hydrology.cee.duke.edu/POLARIS/)
I used mean soil pH in water values collected from 5 to 15 cm deep from the surface, with a resolution of 1 arcsec (~30 meters) 

Nathaniel W. Chaney, Eric F. Wood, Alexander B. McBratney, Jonathan W. Hempel, Travis W. Nauman, Colby W. Brungard, Nathan P. Odgers,
POLARIS: A 30-meter probabilistic soil series map of the contiguous United States, Geoderma, Volume 274, 2016, Pages 54-67, ISSN 0016-7061,https://doi.org/10.1016/j.geoderma.2016.03.025.(https://www.sciencedirect.com/science/article/pii/S0016706116301434)

Insert plot of soil ph

#### Elevation
Elevation was taken from NASA Shuttle Radar Topography Mission (SRTM) provied by the Land Processes Distributed Active Archive Center (LP DAAC). SRTM is a collaboration between NASA, the National Geospatial-Intelligence Agency (NGA), and German and Italian space agencies. It is a near-global digital elevation model (DEM) using radar inferometry. Data was collected on the Space Shuttle Endeavour during its STS-99 mission from 02/11/2000 to 02/22/2000. This is version 3 of SRTM data and has a resolution of 1 arcsec. 

NASA JPL (2013). <i>NASA Shuttle Radar Topography Mission Global 1 arc second</i> [Data set]. NASA EOSDIS Land Processes Distributed Active Archive Center. Accessed 2024-12-15 from https://doi.org/10.5067/MEaSUREs/SRTM/SRTMGL1.003
10

Insert plot of elevation 

#### Precipitation 

Insert plot of precipitation 

### Model
To train a fuzzy logic habitat suitability model:

Research S. nutans, and find out what optimal values are for each variable you are using (e.g. soil pH, slope, and current climatological annual precipitation).
For each digital number in each raster, assign a value from 0 to 1 for how close that grid square is to the optimum range (1=optimal, 0=incompatible).
Combine your layers by multiplying them together. This will give you a single suitability number for each square.
Optionally, you may apply a threshold to make the most suitable areas pop on your map.

https://www.nrcs.usda.gov/plantmaterials/etpmcpg13196.pdf
