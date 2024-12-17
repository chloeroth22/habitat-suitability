# Habitat Suitability 
## Project Description  
In this project, I developed a habitat suitability model for Sorghastrum nutans (commonly known as Indiangrass). This project aims to investigate whether changes in key environmental variables—soil pH, slope, and precipitation—are making conditions less suitable for S. nutans in the South, while making them more favorable in the North, potentially driving the species' northward range expansion. To explore this potential range shift, I compared two grasslands from the USFS National Grasslands. I selected one grassland from the Northern United States and one from the Southern United States to examine possible regional differences. The Northern grassland chosen was the Sheyenne National Grassland in North Dakota. and Southern grassland selected was the Lyndon B. Johnson National Grassland in Texas. 

Soil pH is from from POLARIS, elevation data is from NASA Shuttle Radar Topography Mission (SRTM), and precipitation data is from Bejing Normal University's (BNU-ESM) "Monthly aggregation of downscaled daily meteorological data of Monthly Precipitation Amount". I chose to train a fuzzy logic habitat suitability model. For S. nutans I researched what the optimal values are for each variable (soil pH, slope, and historic annual precipitation) and for each digital number in each raster, a value from 0 to 1 was assigned for how close that grid square is to the optimum range (1=optimal, 0=incompatible).

## 

## Installation
Use the [enviornment earth-analytics-python.yml](earth-analytics-habitat.py)


