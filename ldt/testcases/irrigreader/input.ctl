netcdf ^lis_input.d01 {
dimensions:
	east_west = 381 ;
	north_south = 281 ;
	east_west_b = 385 ;
	north_south_b = 285 ;
	time = 1 ;
	sfctypes = 20 ;
	croptypes = 1 ;
variables:
	float time(time) ;
	float DOMAINMASK(north_south, east_west) ;
		DOMAINMASK:standard_name = "DOMAINMASK" ;
		DOMAINMASK:units = "" ;
		DOMAINMASK:scale_factor = 1.f ;
		DOMAINMASK:add_offset = 0.f ;
		DOMAINMASK:missing_value = -9999.f ;
		DOMAINMASK:vmin = 0.f ;
		DOMAINMASK:vmax = 0.f ;
		DOMAINMASK:num_bins = 1 ;
	float LANDMASK(north_south, east_west) ;
		LANDMASK:standard_name = "LANDMASK" ;
		LANDMASK:units = "" ;
		LANDMASK:scale_factor = 1.f ;
		LANDMASK:add_offset = 0.f ;
		LANDMASK:missing_value = -9999.f ;
		LANDMASK:vmin = 0.f ;
		LANDMASK:vmax = 0.f ;
		LANDMASK:num_bins = 1 ;
	float SURFACETYPE(sfctypes, north_south, east_west) ;
		SURFACETYPE:standard_name = "Surface type" ;
		SURFACETYPE:units = "-" ;
		SURFACETYPE:scale_factor = 1.f ;
		SURFACETYPE:add_offset = 0.f ;
		SURFACETYPE:missing_value = -9999.f ;
		SURFACETYPE:vmin = 0.f ;
		SURFACETYPE:vmax = 0.f ;
		SURFACETYPE:num_bins = 20 ;
	float LANDCOVER(sfctypes, north_south, east_west) ;
		LANDCOVER:standard_name = "MODIS-IGBP (NCEP-modified) landcover map" ;
		LANDCOVER:units = "" ;
		LANDCOVER:scale_factor = 1.f ;
		LANDCOVER:add_offset = 0.f ;
		LANDCOVER:missing_value = -9999.f ;
		LANDCOVER:vmin = 0.f ;
		LANDCOVER:vmax = 0.f ;
		LANDCOVER:num_bins = 20 ;
	float CROPTYPE(north_south, east_west) ;
		CROPTYPE:standard_name = "CROPTYPE" ;
		CROPTYPE:units = "" ;
		CROPTYPE:scale_factor = 1.f ;
		CROPTYPE:add_offset = 0.f ;
		CROPTYPE:missing_value = -9999.f ;
		CROPTYPE:vmin = 0.f ;
		CROPTYPE:vmax = 0.f ;
		CROPTYPE:num_bins = 1 ;
	float IRRIGFRAC(north_south, east_west) ;
		IRRIGFRAC:standard_name = "User Derived Irrig gridcell fraction" ;
		IRRIGFRAC:units = "-" ;
		IRRIGFRAC:scale_factor = 1.f ;
		IRRIGFRAC:add_offset = 0.f ;
		IRRIGFRAC:missing_value = -9999.f ;
		IRRIGFRAC:vmin = 0.f ;
		IRRIGFRAC:vmax = 0.f ;
		IRRIGFRAC:num_bins = 1 ;
	float lat(north_south, east_west) ;
		lat:standard_name = "latitude" ;
		lat:units = "degrees_north" ;
		lat:scale_factor = 1.f ;
		lat:add_offset = 0.f ;
		lat:missing_value = -9999.f ;
		lat:vmin = 0.f ;
		lat:vmax = 0.f ;
	float lon(north_south, east_west) ;
		lon:standard_name = "longitude" ;
		lon:units = "degrees_east" ;
		lon:scale_factor = 1.f ;
		lon:add_offset = 0.f ;
		lon:missing_value = -9999.f ;
		lon:vmin = 0.f ;
		lon:vmax = 0.f ;
	float lat_b(north_south_b, east_west_b) ;
		lat_b:standard_name = "latitude_b" ;
		lat_b:units = "degrees_north" ;
		lat_b:scale_factor = 1.f ;
		lat_b:add_offset = 0.f ;
		lat_b:missing_value = -9999.f ;
		lat_b:vmin = 0.f ;
		lat_b:vmax = 0.f ;
	float lon_b(north_south_b, east_west_b) ;
		lon_b:standard_name = "longitude_b" ;
		lon_b:units = "degrees_east" ;
		lon_b:scale_factor = 1.f ;
		lon_b:add_offset = 0.f ;
		lon_b:missing_value = -9999.f ;
		lon_b:vmin = 0.f ;
		lon_b:vmax = 0.f ;

// global attributes:
		:MAP_PROJECTION = "EQUIDISTANT CYLINDRICAL" ;
		:SOUTH_WEST_CORNER_LAT = 20.975f ;
		:SOUTH_WEST_CORNER_LON = 20.975f ;
		:DX = 0.05f ;
		:DY = 0.05f ;
		:INC_WATER_PTS = "false" ;
		:LANDCOVER_SCHEME = "IGBPNCEP" ;
		:BARESOILCLASS = 16 ;
		:URBANCLASS = 13 ;
		:SNOWCLASS = 15 ;
		:WATERCLASS = 17 ;
		:WETLANDCLASS = 11 ;
		:GLACIERCLASS = 15 ;
		:CROPCLASS = 12 ;
		:NUMVEGTYPES = 17 ;
		:LANDMASK_SOURCE = "MODIS_Native" ;
		:SFCMODELS = "none" ;
		:CROPCLASS_SCHEME = "CROPMAP" ;
		:CROPCLASS_NUMBER = 19 ;
		:SOILTEXT_SCHEME = "Soil texture not selected" ;
		:title = "Land Data Toolkit (LDT) parameter-processed output" ;
		:institution = "NASA GSFC Hydrological Sciences Laboratory" ;
		:history = "created on date: 2020-12-28T15:12:09.345" ;
		:references = "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007" ;
		:comment = "website: http://lis.gsfc.nasa.gov/" ;
}
