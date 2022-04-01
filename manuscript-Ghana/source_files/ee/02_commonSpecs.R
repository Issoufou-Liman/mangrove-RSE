## Common scale for aggregation
AgScale = 30;

## Temporal scale for composite
nDays = 30;

## Lags days dependency
lagDays = nDays + 1;
## print(lagDays)

## The number of cycle per year
harmonics = 2 

## load all functions

## Should local files with the same name be overwritten?
overwrite_file = TRUE

## region of interest
checkPointGhana = ee$Geometry$Point(list(-1.1993780600423731,5.148891150785175));

roiGhana = ee$Geometry$Polygon(
  list(list(list(-3.2795998579167573,5.092430973484513),
            list(-2.0985695844792573,4.71479062410626),
            list(1.1808493608332427,5.852494955748432),
            list(1.1863425248957427,9.277095344549588),
            list(-3.2741066938542573,9.266252547435618),
            list(-3.2795998579167573,5.092430973484513))), NULL, FALSE);

lapply(list('.shp', '.shx', 'prj', 'dbf', 'cpg'), function(i){
  ee_to_drive_to_local (
    ee_object = ee$FeatureCollection(list(ee$Feature(roiGhana))),
    drive_description = '~/Ghana_Manuscript/roiGhana', 
    drive_folder = 'Ghana_Manuscript', 
    region = roiGhana, 
    scale = AgScale,
    maxPixels = 1e13, 
    local_folder  = 'ee_output/roiGhana',
    overwrite = overwrite_file,
    ee_type = 'table',
    file_extension = i
  )
})
## Define a palette for the 17 land cover classes.
myPalette = c(
  'aec3d4', # 1
  '152106', # 2
  '369b47', # 4
  '30eb5b', # 5
  '225129', # 3
  '6a2325', # 6
  'c3aa69', # 7
  '91af40', # 8
  '111149', # 9
  'cdb33b', # 10
  'cc0013', # 11
  '33280d', # 12
  'd7cdcc', # 13
  'ffc0cb', # 14
  'b76031', # 15
  'f7e084', # 16
  'd9903d'  # 17
)

myPalette <- paste0('#', myPalette)

className <- c(
  'Water bodies',
  'Closed forests',
  'Open forests',
  'Woody savannas',
  'Savannas',
  'Closed shrublands',
  'Opened shrublands',
  'Grasslands',
  'Wetlands',
  'Croplands',
  'Built-up areas',
  'Mosaics',
  'Barren',
  'Mangroves',
  'Salt mines',
  'Tree plantations',
  'Riperian vegetation'
)
# maxpixels = 280597463
# maxpixels =   5000000
maxpixels =   5100000

## Change direction ####

dpi = 300

magnitudeClassName = c(
  'Increase',
  'Decrease',
  'Stable',
  'water'
)
myMagnitudePalette <- c(
  'green', # Icrease 
  'red',   # Decrease
  'blue',   # Stable,
  "#D8F1FF" # Water
)

gplot_data <- function(x, maxpixels=ncell(x),...){
  if (maxpixels < ncell(x)){
    x <- sampleRegular(x, size = maxpixels, asRaster=TRUE)
  }
  coords <- xyFromCell(x, seq_len(ncell(x)))
  dat <- stack(as.data.frame(getValues(x)))
  names(dat) <- c('value', 'variable')
  
  dat <- cbind(coords, dat)
}

focal_stack <- function(x){
  n_layers <- seq(nlayers(x))
  r_list <- lapply(X = n_layers, 
                   # function(i) focal(x[[i]], w=matrix(1/9, nc=3, nr=3))
                   function(i) focal(x[[i]], w=matrix(1,3,3), fun=function(j) modal(j, na.rm = FALSE), ties = 'NA', pad=TRUE)
  )
  stack(r_list)
}

