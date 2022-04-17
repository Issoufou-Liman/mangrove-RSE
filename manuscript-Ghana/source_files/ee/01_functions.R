ee_task_monitoring <- function (task, task_time = 5, eeTaskList = FALSE, quiet = FALSE)
{
  if (missing(task)) {
    all_task <- ee_utils_py_to_r(ee$batch$Task$list())
    task <- all_task[[1]]
  }
  if (is.character(task)) {
    all_task <- ee_utils_py_to_r(ee$batch$Task$list())
    id_tasks <- lapply(all_task, function(task) task[["id"]]) %>%
      unlist()
    if (any(id_tasks %in% task)) {
      task <- all_task[[which(id_tasks %in% task)]]
    }
    else {
      stop("Undefined Task ID entered")
    }
  }
  if (eeTaskList) {
    if (!quiet) {
      cat("EETaskList:\n")
    }
    task_list <- mapply(function(x) {
      sprintf("<Task %s: %s (%s)>", x$task_type, x$config,
              x$state)
    }, ee$batch$Task$list())
    if (!quiet) {
      cat("", paste0(task_list, "\n"))
      cat("\n")
    }
  }
  counter <- 0
  while (ee_utils_py_to_r(ee$batch$Task$active(task)) & task[["state"]] !=
         "CANCEL_REQUESTED") {
    if (!quiet) {
      cat(sprintf("Polling for task <id: %s, time: %ds>.\n",
                  task[["id"]], counter))
    }
    counter <- counter + task_time
    Sys.sleep(task_time)
  }
  task_status <- ee_utils_py_to_r(ee$batch$Task$status(task))
  if (!quiet) {
    cat(sprintf("State: %s\n", task_status[["state"]]))
  }
  if (task_status[["state"]] != "COMPLETED") {
    message("ERROR in Earth Engine servers: ", task_status[["error_message"]])
    stop("ee_monitoring was forced to stop before getting results")
  }
  invisible(task)
}

# C:\Users\User\AppData\Local\r-miniconda\envs\rgee\python.exe
# Landsat 8 sanity 
#' Function to mask clouds based on the pixel_qa band of Landsat 8 SR data.
#' 
#' @param {ee.Image} image input Landsat 8 SR image
#' @return {ee.Image} cloudmasked Landsat 8 image
maskL8sr = function (image) {
  # # Bits 3 and 5 are cloud shadow and cloud, respectively.
  cloudShadowBitMask = bitwShiftL (1, 3);
  cloudsBitMask = bitwShiftL (1, 5);
  # # Get the pixel QA band.
  qa = image$select('pixel_qa');
  # Both flags should be set to zero, indicating clear conditions.
  mask = qa$bitwiseAnd(cloudShadowBitMask)$eq(0)$And(qa$bitwiseAnd(cloudsBitMask)$eq(0));
  out <- image$updateMask(mask);
};

# Landsat 7 sanity
#' Function to mask clouds based on the pixel_qa band of Landsat SR data.
#' 
#' @param {ee.Image} image Input Landsat SR image
#' @return {ee.Image} Cloudmasked Landsat image
cloudMaskL457 = function(image) {
  qa = image$select('pixel_qa');
  # If the cloud bit (5) is set and the cloud confidence (7) is high
  # or the cloud shadow bit is set (3), then it's a bad pixel.
  
  cloud = qa$bitwiseAnd(bitwShiftL(1, 5))$And(qa$bitwiseAnd(bitwShiftL(1, 7)))$Or(qa$bitwiseAnd(bitwShiftL(1,3)));
  # Remove edge pixels that don't occur in all bands
  mask2 = image$mask()$reduce(ee$Reducer$min());
  image$updateMask(cloud$Not())$updateMask(mask2);
};

#' Collect images from LANDSAT 4,5,7 and 8 based spatial and temporal filter 
#' @param  {String} collection Either 'LSR' (default) or 'TAO' indicating the collection of interest
#' @param  {FeatureCollection} roi The area of interset
#' @param  {Integer} startYear The beginning of the collection (e.g. 2000)
#' @param  {Integer} endYear The end of the collection
#' @param  {Integer} months Integer between 1 to 12 indicating the months of interest
#' @param  {Integer} nDays Integer between indicating the number of days to consider for composite
#' @param  {Logical} includeL7 Whether to include LANDSAT 7 data in the collection. default to false
#' @return {ImageCollection}      An image collection
collectLansatImages = function  (collection = 'LSR', 
                                 roi = ee$Geometry$Point(list(-1.1993780600423731,5.148891150785175)),
                                 startYear = 1970, 
                                 endYear = 2020,
                                 months = ee$List$sequence(1, 12),
                                 nDays = 16,
                                 includeL7 = FALSE
){
  
  startDate = ee$Date$fromYMD(startYear, 1, 1);
  endDate = ee$Date$fromYMD(endYear, 12, 31);
  
  if (collection == 'LSR'){
    # Assign a common name to the sensor-specific bands.
    l8_bands = list('B2',   'B3',    'B4',  'B5',  'B6',    'B7',    'B10'); #Landsat 8
    l7_bands = list('B1',   'B2',    'B3',  'B4',  'B5',    'B7',    'B6'); #Landsat 7
    l5_bands = list('B1',   'B2',    'B3',  'B4',  'B5',    'B7',    'B6'); #Landsat 4 and 5
    band_col = list('blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'temp');
    
    # Landsat 8 sanity 
    #' Function to mask clouds based on the pixel_qa band of Landsat 8 SR data.
    #' 
    #' @param {ee.Image} image input Landsat 8 SR image
    #' @return {ee.Image} cloudmasked Landsat 8 image
    maskL8sr = function (image) {
      # # Bits 3 and 5 are cloud shadow and cloud, respectively.
      cloudShadowBitMask = bitwShiftL (1, 3);
      cloudsBitMask = bitwShiftL (1, 5);
      # # Get the pixel QA band.
      qa = image$select('pixel_qa');
      # Both flags should be set to zero, indicating clear conditions.
      mask = qa$bitwiseAnd(cloudShadowBitMask)$eq(0)$And(qa$bitwiseAnd(cloudsBitMask)$eq(0));
      out <- image$updateMask(mask);
    };
    
    # Landsat 7 sanity
    #' Function to mask clouds based on the pixel_qa band of Landsat SR data.
    #' 
    #' @param {ee.Image} image Input Landsat SR image
    #' @return {ee.Image} Cloudmasked Landsat image
    cloudMaskL457 = function(image) {
      qa = image$select('pixel_qa');
      # If the cloud bit (5) is set and the cloud confidence (7) is high
      # or the cloud shadow bit is set (3), then it's a bad pixel.
      
      cloud = qa$bitwiseAnd(bitwShiftL(1, 5))$And(qa$bitwiseAnd(bitwShiftL(1, 7)))$Or(qa$bitwiseAnd(bitwShiftL(1,3)));
      # Remove edge pixels that don't occur in all bands
      mask2 = image$mask()$reduce(ee$Reducer$min());
      image$updateMask(cloud$Not())$updateMask(mask2);
    };
    
    img_col = ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$
      merge(ee$ImageCollection('LANDSAT/LC08/C01/T2_SR'))$
      filterDate(startDate, endDate)$
      filterBounds(roi)$
      map(ee_utils_pyfunc(maskL8sr))$
      select(l8_bands, band_col);
    
    if(includeL7){
      img_col = ee$ImageCollection(img_col$merge(
        ee$ImageCollection('LANDSAT/LE07/C01/T1_SR')$
          merge(ee$ImageCollection('LANDSAT/LE07/C01/T2_SR'))$
          filterDate(startDate, endDate)$
          filterBounds(roi)$
          map(ee_utils_pyfunc(cloudMaskL457))$
          select(l7_bands, band_col)
      ));
    }
    
    img_col = ee$ImageCollection(img_col$merge(
      ee$ImageCollection('LANDSAT/LT05/C01/T1_SR')$
        merge(ee$ImageCollection('LANDSAT/LT05/C01/T2_SR'))$
        filterDate(startDate, endDate)$
        filterBounds(roi)$
        map(ee_utils_pyfunc(cloudMaskL457))$
        select(l5_bands, band_col)
    )); 
    
    img_col = ee$ImageCollection(img_col$merge(
      ee$ImageCollection('LANDSAT/LT04/C01/T1_SR')$
        merge(ee$ImageCollection('LANDSAT/LT04/C01/T2_SR'))$
        filterDate(startDate, endDate)$
        filterBounds(roi)$
        map(ee_utils_pyfunc(cloudMaskL457))$
        select(l5_bands, band_col)
    ));
    
  } else if (collection == 'TOA'){
    # Function to mask clouds in Landsat 8 imagery.
    maskClouds = function(image) {
      quality = image$select('BQA');
      cloud01 = quality$eq(61440);
      cloud02 = quality$eq(53248);
      cloud03 = quality$eq(28672);
      mask = cloud01$or(cloud02)$or(cloud03)$not();
      return (image$updateMask(mask));
    };
    
    # Assign a common name to the sensor-specific bands.
    l8_bands = list('B2',   'B3',    'B4',  'B5',  'B6',    'B7',    'B10'); #Landsat 8
    l7_bands = list('B1',   'B2',    'B3',  'B4',  'B5',    'B7',    'B6_VCID_1'); #Landsat 7
    l5_bands = list('B1',   'B2',    'B3',  'B4',  'B5',    'B7',    'B6'); #Landsat 4 and 5
    band_col = list('blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'temp');
    
    img_col = ee$ImageCollection('LANDSAT/LC08/C01/T1_TOA')$
      filterDate(startDate, endDate)$
      filterBounds(roi)$
      # map(ee_utils_pyfunc(maskClouds))$
      select(l8_bands, band_col);
    if(includeL7){
      img_col = ee$ImageCollection(img_col$merge(
        ee$ImageCollection('LANDSAT/LE07/C01/T1_TOA')$
          filterDate(startDate, endDate)$
          filterBounds(roi)$
          # map(ee_utils_pyfunc(maskClouds))$
          select(l7_bands, band_col)
      ));
    }            
    
    img_col = ee$ImageCollection(img_col$merge(
      ee$ImageCollection('LANDSAT/LT05/C01/T1_TOA')$
        filterDate(startDate, endDate)$
        filterBounds(roi)$
        # map(ee_utils_pyfunc(maskClouds))$
        select(l5_bands, band_col)
    ));
    
    img_col = ee$ImageCollection(img_col$merge(
      ee$ImageCollection('LANDSAT/LT04/C01/T1_TOA')$
        filterDate(startDate, endDate)$
        filterBounds(roi)$
        # map(ee_utils_pyfunc(maskClouds))$
        select(l5_bands, band_col)
    ));
    
  }
  # limit the collection to some production date for reproducing the same result
  img_col = img_col$
    filterMetadata(
      'LEVEL1_PRODUCTION_DATE',
      'less_than',
      # ee$Date('2021-03-11T00:00:00')$millis() #  '2021-03-11T05:45:16'
      # ee$Date('2021-03-12T19:00:00')$millis()
      # ee$Date('2021-03-12T23:59:59')$millis()
      ee$Date('2021-03-08T21:45:32')$millis() # 210308_2 attempt
    );
  # Making n days composite
  datesList = ee$List(img_col$aggregate_array('system:time_start'))$
    distinct()$
    map(ee_utils_pyfunc(function(dateMillis){
      return (ee$Date(dateMillis));
    }));
  
  # Recomputed startDate and endDate based on the available images in the collection
  startDate = ee$Date(datesList$sort()$get(0)); 
  endDate = ee$Date(datesList$sort()$get(datesList$size()$subtract(2)));
  
  imageDates = ee$List$sequence(0, endDate$difference(startDate,'day')$round(), nDays)$
    map(ee_utils_pyfunc(function(dy) {
      return (startDate$advance(dy,'day'));
    }));
  
  imageOut = imageDates$map(ee_utils_pyfunc(function(dt){
    out = img_col$filterDate(ee$Date(dt), ee$Date(dt)$advance(nDays-1, 'day'));
    out = out$median();
    
    out = ee$Algorithms$If(
      condition = out$bandNames()$size()$neq(0),
      trueCase = ee$Image(out),
      falseCase = ee$Image(list(0,0,0,0,0,0,0))$selfMask()$rename(band_col)
    );
    
    out = ee$Image(out)$clip(roi);
    
    # Compute time in fractional years since the epoch.
    years = ee$Date(dt)$difference(ee$Date('1970-01-01'), 'year');
    
    # comopute the normalized difference spectal index
    GoaNDWI      = out$normalizedDifference(list('nir', 'swir1'))$rename('NDWI');
    NDVI         = out$normalizedDifference(list('nir', 'red'))$rename('NDVI');
    
    # Return the image with the added bands.
    return (
      out$# select('temp')
        addBands(GoaNDWI)$
        addBands(NDVI)$
        
        # Add a time band.
        addBands(ee$Image(years)$rename('t'))$
        float()$
        
        # Add a constant band.
        addBands(ee$Image$constant(1))$
        set('system:time_start', ee$Date(dt)$millis())
    );
  }));
  return (imageOut);
}

lag = function(leftCollection, rightCollection, lagDays) {
  filter = ee$Filter$And(
    ee$Filter$maxDifference(
      difference = 1000 * 60 * 60 * 24 * lagDays,
      leftField = 'system:time_start', 
      rightField = 'system:time_start'
    ), 
    ee$Filter$greaterThan(
      leftField = 'system:time_start', 
      rightField = 'system:time_start'
    ));
  
  return (ee$Join$saveAll(
    matchesKey = 'images',
    measureKey = 'delta_t',
    ordering = 'system:time_start', 
    ascending = FALSE ## Sort reverse chronologically
  )$apply(
    primary = leftCollection, 
    secondary = rightCollection, 
    condition = filter
  ));
}

merge = function(image) {
  # Function to be passed to iterate.
  merger = function(current, previous) {
    return (ee$Image(previous)$addBands(current));
  };
  return (ee$ImageCollection$fromImages(
    image$get('images'))$iterate(merger, image));
}

covariance = function(mergedCollection, band, lagBand) {
  return (
    mergedCollection$select(list(band, lagBand))$map(ee_utils_pyfunc(function(image) {
      return (image$toArray());
    }))$reduce(ee$Reducer$covariance(), 16)
  );
}

correlation = function(vcArrayImage) {
  covariance = ee$Image(vcArrayImage)$arrayGet(list(0, 1));
  sd0 = ee$Image(vcArrayImage)$arrayGet(list(0, 0))$sqrt();
  sd1 = ee$Image(vcArrayImage)$arrayGet(list(1, 1))$sqrt();
  return (covariance$divide(sd0)$divide(sd1)$rename('correlation'));
}


harmonic_fill = function(imgCollection, dependent = imgCollection$first()$bandNames(), harmonics = 1, dateOrigin = '1970-01-01'){
  # Zhu and woodcock min number of observations
  minNumObsMask = imgCollection$count()$gte(15);
  # return (minNumObsMask);
  
  # tmaskBands = tmaskBands || ee$List(list('green', 'swir1')) 
  # dependent = dependent || imgCollection$first()$bandNames();
  # Function to get a sequence of band names for harmonic terms.
  # getNames = function(base, list) {
  #   return (ee$List(list)$map(ee_utils_pyfunc(function(i) {
  #     return (ee$String(base)$cat(ee$String(i)$int()));
  #   })));
  # };
  
  getNames = function(base, list) {
    out <- sapply(list, function(x){
      paste0(base, x)
    })
  };
  
  # Make a list  and band name suffixes out of harmonics .  
  # harmonicFrequencies = ee$List$sequence(1, harmonics);
  harmonicFrequencies = 1:harmonics;
  # Construct lists of names for the harmonic terms.
  cosNames = getNames('cos_', harmonicFrequencies);
  sinNames = getNames('sin_', harmonicFrequencies);
  # Independent variables.
  baseIndependents = ee$List(list('constant', 't'));
  independents = baseIndependents$cat(cosNames)$cat(sinNames);
  
  imgCollection = imgCollection$map(ee_utils_pyfunc(function(img){
    img = ee$Image(img)$updateMask(minNumObsMask);
    # get the image date
    date = ee$Date(img$get('system:time_start'));
    # get the fractional year since dateOrigin
    years = date$difference(ee$Date(dateOrigin), 'year');
    # Convert the fractional years to radian and cast it to ee$Image
    timeRadians = ee$Image(years$multiply(2 * pi));
    # Make an image out of the frequencies.
    frequencies = ee$Image$constant(harmonicFrequencies);
    # Get the cosine terms.
    cosines = timeRadians$multiply(frequencies)$cos()$rename(cosNames);
    # Get the sin terms.
    sines = timeRadians$multiply(frequencies)$sin()$rename(sinNames);
    return (img$addBands(ee$Image(1))$
              addBands(timeRadians$rename('t')$float())$
              addBands(cosines)$
              addBands(sines));
  }));
  # The output of the regression reduction is an array image.
  harmonicTrend = imgCollection$
    select(independents$cat(dependent))$
    reduce(ee$Reducer$robustLinearRegression(independents$length(), dependent$length()));
  
  harmonicTrendCoefficients = harmonicTrend$select('coefficients')
  
  cosCoeffArray = harmonicTrendCoefficients$
    arraySlice(
      axis = 0, 
      start = baseIndependents$size(), 
      end = baseIndependents$size()$add(harmonics),
      step = 1
    );
  
  sinCoeffArray = harmonicTrendCoefficients$
    arraySlice(
      axis = 0, 
      start = baseIndependents$size()$add(harmonics), 
      end =  baseIndependents$size()$add(harmonics)$add(harmonics),
      step = 1
    );
  
  # Compute phase and amplitude.
  phase = sinCoeffArray$atan2(cosCoeffArray)
  
  amplitude = sinCoeffArray$hypot(cosCoeffArray);
  
  fittedBandName = dependent$map(ee_utils_pyfunc(function(dep){
    return (ee$String(dep)$cat('_')$cat('fitted'));
  }));
  filledBandName = dependent$map(ee_utils_pyfunc(function(dep){
    return (ee$String(dep)$cat('_')$cat('filled'));
  }));
  
  # Get and scale rmse by 3 (equation 4 zhu and Woodcock, 2014)
  residualsBandName = dependent$map(ee_utils_pyfunc(function(dep){
    return (ee$String(dep)$cat('_')$cat('residuals'));
  }));
  rmse = harmonicTrend$select('residuals');
  
  fittedHarmonic = imgCollection$map(ee_utils_pyfunc(function(image) {
    # assuming that number of independents = 6 and number of bands = 2,
    # harmonicTrendCoefficients is a 6 x 2 array
    indepArray = image$select(independents)$toArray()$toArray(1) # 6 x 1
    depArray = image$select(dependent)$toArray()$toArray(0)
    # transpose indepArray to  1 x 6: number of column should be 
    # equal to number of rows (6 x 2) in harmonicTrendCoefficients for compatibility in matrix multiplication
    fitted = indepArray$arrayTranspose()$
      matrixMultiply(harmonicTrendCoefficients)$ # 1 x 2
      arrayReduce(
        reducer = ee$Reducer$sum(), 
        axes = list(0)
      )$
      arrayProject(list(1));
    
    zhu_and_woodcock_B2Mask = depArray$arraySlice(axis = 0, start = 1, end = 2, step = 1)$
      subtract(
        fitted$arraySlice(axis = 0, start = 1, end = 2, step = 1)
      )$gt(0.04)$
      Not();
    
    zhu_and_woodcock_B5Mask = depArray$arraySlice(axis = 0, start = 4, end = 5, step = 1)$
      subtract(
        fitted$arraySlice(axis = 0, start = 4, end = 5, step = 1)
      )$lt(-0.04)$
      Not();
    
    zhu_and_woodcock_Mask   = zhu_and_woodcock_B2Mask$Or(zhu_and_woodcock_B5Mask)$arrayFlatten(list(list('zhu_and_woodcock_Mask')));
    depArray = depArray$updateMask(zhu_and_woodcock_Mask);
    
    filled = fitted$where(depArray$mask(), depArray)#$rename(filledBandName);
    
    out = image$select(dependent)$
      addBands(filled$arrayFlatten(list(filledBandName)))$
      addBands(fitted$arrayFlatten(list(fittedBandName)));
    
    # zhu_and_woodcock_B2Mask = out$select('green')$subtract(out$select('green_fitted'))$gt(0.04)$Not();
    # zhu_and_woodcock_B5Mask = out$select('swir1')$subtract(out$select('swir1_fitted'))$lt(-0.04)$Not()
    # zhu_and_woodcock_Mask   = zhu_and_woodcock_B2Mask$Or(zhu_and_woodcock_B5Mask);
    
    out = out$ #$updateMask(zhu_and_woodcock_Mask)
      set('system:time_start', image$get('system:time_start'));
    return (out);
    
  }));
  
  return (fittedHarmonic$
            set('rmse', rmse$arrayProject(list(0))$arrayFlatten(list(residualsBandName)))$
            set('phase', phase)$
            set('amplitude', amplitude))
  ;
}


#'Fit harmonic trend to image collection
#'@param  {ee.ImageCollection} imgCollection The collection of image. Must contained a valide time field
#'@param  {String} dateOrigin the date origin from which to compute the franctional years. Default to the epoch ('1970-01-01')
#'@param  {Integer} harmonics the number of cyles to be considered. Default to 1 cycle
#'@param  {String} the dependent band variable
#'@return {ImageCollection}      A collection of image bands
harmonic_fit = function(imgCollection, dependent = imgCollection$first()$bandNames(), harmonics = 1, lagDays, dateOrigin = '1970-01-01'){
  # Function to get a sequence of band names for harmonic terms.
  # getNames = function(base, list) {
  #   return (ee$List(list)$map(ee_utils_pyfunc(function(i) { 
  #     return (ee$String(base)$cat(ee$Number(i)$int()));
  #   })));
  # };
  # 
  getNames = function(base, list) {
    out <- sapply(list, function(x){
      paste0(base, x)
    })
  };
  # Make a list  and band name suffixes out of harmonics .  
  # harmonicFrequencies = ee$List$sequence(1, harmonics);
  harmonicFrequencies = 1:harmonics;
  # Construct lists of names for the harmonic terms.
  cosNames = getNames('cos_', harmonicFrequencies);
  sinNames = getNames('sin_', harmonicFrequencies);
  # Independent variables.
  baseIndependents = ee$List(list('constant', 't'));
  independents = baseIndependents$cat(cosNames)$cat(sinNames);
  
  imgCollection = imgCollection$map(ee_utils_pyfunc(function(img){
    img = ee$Image(img);
    # get the image date
    date = ee$Date(img$get('system:time_start'));
    # get the fractional year since dateOrigin
    years = date$difference(ee$Date(dateOrigin), 'year');
    # Convert the fractional years to radian and cast it to ee$Image
    timeRadians = ee$Image(years$multiply(2 * pi));
    
    # Make an image out of the frequencies.
    frequencies = ee$Image$constant(harmonicFrequencies);
    
    # Get the cosine terms.
    cosines = timeRadians$multiply(frequencies)$cos()$rename(cosNames);
    # Get the sin terms.
    sines = timeRadians$multiply(frequencies)$sin()$rename(sinNames);
    
    return (
      img$addBands(ee$Image(1))$
        addBands(timeRadians$rename('t')$float())$
        addBands(cosines)$
        addBands(sines)
    );
  }));
  # The output of the regression reduction is an array image.
  harmonicTrend = imgCollection$
    select(independents$cat(dependent))$
    reduce(ee$Reducer$linearRegression(independents$length(), dependent$length()));
  
  harmonicTrendCoefficients = harmonicTrend$select('coefficients')
  
  cosCoeffArray = harmonicTrendCoefficients$
    arraySlice(
      axis = 0, 
      start = baseIndependents$size(), 
      end = baseIndependents$size()$add(harmonics),
      step = 1
    );
  
  sinCoeffArray = harmonicTrendCoefficients$
    arraySlice(
      axis = 0, 
      start = baseIndependents$size()$add(harmonics), 
      end =  baseIndependents$size()$add(harmonics)$add(harmonics),
      step = 1
    );
  
  # Compute phase and amplitude.
  phase = sinCoeffArray$atan2(cosCoeffArray)
  
  amplitude = sinCoeffArray$hypot(cosCoeffArray);
  
  fittedBandName = dependent$map(ee_utils_pyfunc(function(dep){
    return (ee$String(dep)$cat('_')$cat('fitted'));
  }));
  filledBandName = dependent$map(ee_utils_pyfunc(function(dep){
    return (ee$String(dep)$cat('_')$cat('filled'));
  }));
  
  # Get and scale rmse by 3 (equation 4 zhu and Woodcock, 2014)
  residualsBandName = dependent$map(ee_utils_pyfunc(function(dep){
    return (ee$String(dep)$cat('_')$cat('residuals'));
  }));
  rmse = harmonicTrend$select('residuals')#$arrayProject(list(0))$arrayFlatten(list(residualsBandName));
  # rmse = rmse$multiply(ee$Number(1)$divide(3))
  
  fittedHarmonic = imgCollection$map(ee_utils_pyfunc(function(image) {
    # assuming that number of independents = 6 and number of bands = 2,
    # harmonicTrendCoefficients is a 6 x 2 array
    indepArray = image$select(independents)$toArray()$toArray(1) # 6 x 1
    depArray = image$select(dependent)$toArray()$toArray(0)
    # transpose indepArray to  1 x 6: number of column should be 
    # equal to number of rows (6 x 2) in harmonicTrendCoefficients for compatibility in matrix multiplication
    fitted = indepArray$arrayTranspose()$
      matrixMultiply(harmonicTrendCoefficients)$ # 1 x 2
      arrayReduce(
        reducer = ee$Reducer$sum(), 
        axes = list(0)
      )$
      arrayProject(list(1));
    
    filled = fitted$where(depArray$mask(), depArray)#$rename(filledBandName);
    return (
      image$select(dependent$cat(baseIndependents))$
        addBands(filled$arrayFlatten(list(filledBandName)))$
        addBands(fitted$arrayFlatten(list(fittedBandName)))$
        set('system:time_start', image$get('system:time_start'))
    );
    
  }));
  
  fittedHarmonic = fittedHarmonic$
    # $set('rmse', rmse$arrayProject(list(0))$arrayFlatten(list(residualsBandName)))
    set('phase', phase)$
    set('amplitude', amplitude)
  ;
  
  # linear trend 
  # List of the independent variable names
  indepLinear = ee$List(list('constant', 't'));
  linearFitBandNames = dependent$map(ee_utils_pyfunc(function(dep){
    return (ee$String(dep)$cat('_')$cat('linearTrend'));
  }));
  
  # Compute the linear trend. The computation returns the residuals an coefficients (2 x 1) as bands.
  trend = fittedHarmonic$select(indepLinear$cat(dependent))$
    reduce(ee$Reducer$linearRegression(indepLinear$length(), dependent$length()));
  
  # Flatten the coefficients into a 2-band image
  coefficients = trend$select('coefficients')
  #   # $arrayProject(list(0))
  #   $arrayFlatten(list(indepLinear, dependent));
  # return coefficients
  # Compute the trend of the series.
  fittedLinear = fittedHarmonic$map(ee_utils_pyfunc(function(image) {
    # return image
    # $addBands(
    #   image$select(dependent)
    #   $subtract(image$select(indepLinear)$multiply(coefficients)$reduce('sum'))$rename('LinearyDetrented')
    #   $copyProperties(image, list('system:time_start'))
    #   );
    indepArray = image$select(indepLinear)$toArray()$toArray(1) # 6 x 1
    depArray = image$select(dependent)$toArray()$toArray(0)
    # transpose indepArray to  1 x 6: number of column should be 
    # equal to number of rows (6 x 2) in harmonicTrendCoefficients for compatibility in matrix multiplication
    fitted = indepArray$arrayTranspose()$
      matrixMultiply(coefficients)$ # 1 x 2
      arrayReduce(
        reducer = ee$Reducer$sum(), 
        axes = list(0)
      )$
      arrayProject(list(1))$
      arrayFlatten(list(linearFitBandNames));
    
    return (fitted$copyProperties(image, list('system:time_start')))
  }));
  return (
    fittedHarmonic$combine(fittedLinear)$
      set('LinearCoefficients', coefficients$arrayProject(list(0))$arrayFlatten(list(dependent)))$
      # $set('linearTrend', fittedLinear)
      set('rmse', rmse$arrayProject(list(0))$arrayFlatten(list(residualsBandName)))$
      
      set('phase', phase$arrayProject(list(0))$arrayFlatten(list(dependent)))$
      set('amplitude', amplitude$arrayProject(list(0))$arrayFlatten(list(dependent)))
    # $set('covariance', covariance17);
  )
}

#' Classify an image using unsupervised classification
#' @param  {Image} image The image to classify
#' @param  {Integer} nClass  The number of class to consider
#' @param  {Integer} scale   The resolution in meters to consider when sampling the training data
#' @param  {FeatureCollection} roi The region of interest
#' @return {ImageCollection}      A classifed images
unsuperClassifyThis = function (image, nClass = 15, scale = 30, roi, tileScale = 1, numPixels = 500){
  # Make the training dataset.
  training = image$sample(
    region = roi,
    scale = scale,
    numPixels = numPixels,
    tileScale = tileScale
  );
  # Instantiate the clusterer and train it.
  clusterer = ee$Clusterer$wekaKMeans(nClass)$train(training);
  # Cluster the input using the trained clusterer.
  return (image$cluster(clusterer));
}

#' Classify an image using training data
#' @param  {Image} image The image to classify
#' @param  {FeatureCollection} training The training data based on which to classify the image
#' @param  {String} propertyName The string name of the property in the training data
#' @param  {String} outputMode Either 'CLASSIFICATION' (default) indicating to return the classes or 'PROBABILITY' indicating to return class probabilities
#' @param  {Integer} numberOfTrees The number (default to 10) of decision trees to create.
#' @return {ImageCollection}      A classifed image
superClassifyThis = function  (image, training, propertyName = 'landcover', outputMode = 'CLASSIFICATION', numberOfTrees = 10){
  # Get a list of bands for supervized classification
  bands = image$bandNames();
  if (training$geometry()$length() != 0){
    if (outputMode == 'CLASSIFICATION'){
      # Train the classifier for classification
      classifier = ee$Classifier$smileRandomForest(numberOfTrees)$setOutputMode('CLASSIFICATION')$train(
        features = training,
        classProperty = propertyName,
        inputProperties = bands
      );
    } else if (outputMode == 'PROBABILITY'){
      # Train the classifier for classification
      classifier = ee$Classifier$smileRandomForest(numberOfTrees)$setOutputMode('PROBABILITY')$train(
        features = training,
        classProperty = propertyName,
        inputProperties = bands
      );
    }
  }
  # Run the classification
  classified = image$select(bands)$classify(classifier);
  return (classified);
}


#' Assess the accuracy of a classifier based on validation data
#' @param  {Image} image The image based on which the classifier was trained
#' @param  {classifier} classifier The optinal classifier to be used for classifying the image
#' @param  {FeatureCollection} validation The validation data based on which to validate the classification
#' @param  {String} actual The name of the property containing the actual value.
#' @param  {String} predicted The name of the property containing the predicted value.
#' @return {FeatureCollection}      A FeatureCollection containing the Resubstitution error matrix, the overall training accuracy, the validation error matrix, the overall validation accuracy
validateClassification = function (image, training, validation, actual = 'landcover', predicted = 'classification', numberOfTrees){
  bands = image$bandNames();
  # Train the classifier for classification
  classifier = ee$Classifier$smileRandomForest(numberOfTrees)$setOutputMode('CLASSIFICATION')$train(
    features = training,
    classProperty = actual,
    inputProperties = bands
  );
  #  Get a confusion matrix representing resubstitution accuracy.
  trainAccuracy = classifier$confusionMatrix();
  # Classify the validation data.
  validated = validation$classify(classifier);
  # Get a confusion matrix representing expected accuracy.
  testAccuracy = validated$errorMatrix(actual, predicted);
  accuracy = ee$Feature(
    NULL,
    c(
      'Resubstitution error matrix' = trainAccuracy,
      'Validation error matrix' = testAccuracy,
      'Producer accuracy' = testAccuracy$producersAccuracy()$toList()$flatten(),
      'Consumers accuracy' = testAccuracy$consumersAccuracy()$toList()$flatten(),
      'Training overall accuracy' = trainAccuracy$accuracy(),
      'Validation overall accuracy' = testAccuracy$accuracy(),
      'Kappa'  = testAccuracy$kappa()
    )
  );
  return (accuracy);
}

#' Collect images from LANDSAT 4,5,7 and 8 based spatial and temporal filter 
#' 
#' @param  {String} collection Either 'LSR' (default) or 'TAO' indicating the collection of interest
#' @param  {FeatureCollection} roi The area of interset
#' @param  {Integer} startYear The beginning of the collection (e.g. 2000)
#' @param  {Integer} endYear The end of the collection
#' @param  {Integer} months Integer between 1 to 12 indicating the months of interest
#' @param  {Logical} includeL7 Whether to include LANDSAT 7 data in the collection. default to false
#' @return {ImageCollection}      A collection of image bands
collectLansatImages_0 <- function(collection = 'LSR', roi, startYear, endYear, months = ee.List.sequence(1, 12), nDays = 16, includeL7 = FALSE){
  startDate = ee$Date$fromYMD(startYear, 1, 1);
  endDate = ee$Date$fromYMD(endYear, 12, 31);
  if (collection == 'LSR'){
    ## Assign a common name to the sensor-specific bands.
    l8_bands = list('B2',   'B3',    'B4',  'B5',  'B6',    'B7',    'B10'); ##Landsat 8
    l7_bands = list('B1',   'B2',    'B3',  'B4',  'B5',    'B7',    'B6'); ##Landsat 7
    l5_bands = list('B1',   'B2',    'B3',  'B4',  'B5',    'B7',    'B6'); ##Landsat 4 and 5
    modis_bands = list('sur_refl_b03', 'sur_refl_b04', 'sur_refl_b01', 'sur_refl_b02', 'sur_refl_b06', 'sur_refl_b07', 'LST_Day_1km')
    band_col = list('blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'temp');
    
    img_col = ee$ImageCollection('LANDSAT/LC08/C01/T1_SR')$
      merge(ee$ImageCollection('LANDSAT/LC08/C01/T2_SR'))$
      filterDate(startDate, endDate)$
      filterBounds(roi)$
      map(ee_utils_pyfunc(maskL8sr))$
      select(l8_bands, band_col);
    
    if(includeL7){
      img_col = ee$ImageCollection(img_col$merge(
        ee$ImageCollection('LANDSAT/LE07/C01/T1_SR')$
          merge(ee$ImageCollection('LANDSAT/LE07/C01/T2_SR'))$
          filterDate(startDate, endDate)$
          filterBounds(roi)$
          map(ee_utils_pyfunc(cloudMaskL457))$
          select(l7_bands, band_col)
      ));
    }
    
    img_col = ee$ImageCollection(img_col$merge(
      ee$ImageCollection('LANDSAT/LT05/C01/T1_SR')$
        merge(ee$ImageCollection('LANDSAT/LT05/C01/T2_SR'))$
        filterDate(startDate, endDate)$
        filterBounds(roi)$
        map(ee_utils_pyfunc(cloudMaskL457))$
        select(l5_bands, band_col)
    )); 
    
    img_col = ee$ImageCollection(img_col$merge(
      ee$ImageCollection('LANDSAT/LT04/C01/T1_SR')$
        merge(ee$ImageCollection('LANDSAT/LT04/C01/T2_SR'))$
        filterDate(startDate, endDate)$
        filterBounds(roi)$
        map(ee_utils_pyfunc(cloudMaskL457))$
        select(l5_bands, band_col)
    ));
  } 
  return (img_col);
  
}

getAccuracyElement = function(accuracy){
  consumersAccuracy = ee$List(accuracy$get('Consumers accuracy'));
  producerAccuracy = ee$List(accuracy$get('Producer accuracy'));
  
  Kappa = ee$List(accuracy$get('Kappa'));
  overall = ee$List(accuracy$get('Validation overall accuracy'));
  
  klasses  = ee$List$sequence(0, consumersAccuracy$size()$subtract(1))
  return (
    ee$FeatureCollection(
      klasses$map(ee_utils_pyfunc(function(element){
        return (
          ee$Feature(NULL)$
            set('Consumers accuracy', consumersAccuracy$get(element))$
            set('Producer accuracy', producerAccuracy$get(element))$
            set('Kappa', Kappa)$
            set('Overall accuracy', overall)
        );
      }))
    )
  );
}

classAcreages = function(landUse, region, scale, tileScale){
  pixelArea = ee$Image$pixelArea()$
    reproject(
      crs ='EPSG:4326', 
      scale = scale
    )$
    divide(10000); # in hectares
  
  landUse = landUse$mask()$multiply(pixelArea)$addBands(landUse)$
    reduceRegion(
      reducer = ee$Reducer$sum()$group(
        groupField = 1,
        groupName = 'class code'
      ),
      geometry = region,
      scale = scale,
      maxPixels = 1e13,
      tileScale = tileScale
    );
  
  return (
    ee$List(landUse$get('groups'))$
      map(ee_utils_pyfunc(function(obj){
        return (
          ee$Feature(NULL, obj)
          # $set(id = 'landUse');
        )
      }))
  );
}

classDynamics = function(past, present, region, scale, tileScale, targetClass){
  targetClass = ee$Number(targetClass)
  pixelArea = ee$Image$pixelArea()$
    reproject(
      crs ='EPSG:4326', 
      scale = scale
    )$
    divide(10000); # in hectares
  
  # class gain
  gain = past$neq(targetClass)$And(present$eq(targetClass));
  gain = gain$multiply(pixelArea)$addBands(past);
  gain_sums = gain$reduceRegion(
    reducer = ee$Reducer$sum()$group(
      groupField = 1,
      groupName = 'code'
    ),
    geometry = region,
    scale = scale,
    maxPixels = 1e13,
    tileScale = tileScale
  );
  # class loss
  loss = past$eq(targetClass)$And(present$neq(targetClass));
  loss = loss$multiply(pixelArea)$addBands(present);
  loss_sums = loss$reduceRegion(
    reducer = ee$Reducer$sum()$group(
      groupField = 1,
      groupName = 'code'
    ),
    geometry = region,
    scale = scale,
    maxPixels = 1e13,
    tileScale = tileScale
  );
  
  elts = ee$List(gain_sums$get('groups'));
  return (
    ee$List$sequence(0, elts$size()$subtract(1))$
      map(ee_utils_pyfunc(function(id){
        gain = ee$Dictionary(ee$List(gain_sums$get('groups'))$get(id));
        loss = ee$Dictionary(ee$List(loss_sums$get('groups'))$get(id));
        # rename the dictionary to avoid overwriting data as loss and gain have the same keys. 
        renameDictionary = function(dico, prefix){
          oldNames = dico$keys();
          newNames = oldNames$map(ee_utils_pyfunc(function(key){
            return (ee$String(prefix)$cat('_')$cat(key));
          }));
          return (dico$rename(oldNames, newNames));
        };
        
        return (
          ee$Feature(NULL)$
            set(renameDictionary(gain, 'gain'))$
            set(renameDictionary(loss, 'loss'))
        );
      }))
  );
}

classChangeMagnitude = function(past, present, targetClass){
  targetClass = ee$Number(targetClass);
  
  increase = ee$Image(0)$toFloat()$where(
    past$neq(targetClass)$
      And(
        present$eq(targetClass)
      ), 1
  );
  
  decrease =  ee$Image(0)$toFloat()$where(
    past$eq(targetClass)$
      And(
        present$neq(targetClass)
      ), 2
  );
  
  stable =  ee$Image(0)$toFloat()$where(
    past$eq(targetClass)$
      And(
        present$eq(targetClass)
      ), 3
  );  
  
  magnitude = ee$ImageCollection(list(increase, decrease, stable))$
    sum()$
    selfMask();
  
  return (magnitude);    
}


asset_exists <- function(assetid, folder = NULL){
  if(is.null(folder)){
    folder <- ee_get_assethome()
  }
  # asset_path <- paste0(folder, "/", assetid)
  cmd <- paste('earthengine ls', asset_path)
  out <- system(cmd, intern = TRUE)
  out <- sapply(out, function(i){
    fil <- unlist(strsplit(i, '/'))
    fil <- fil[length(fil)]
  })
  assetid %in% out
}

drive_file_exists <- function(target_file){
  short_descrip = tail(unlist(strsplit(target_file, '/')), n =1)
  base_path = gsub(pattern = paste0('/', short_descrip), replacement = '', x = target_file)
  
  file_list <- drive_ls(path = base_path, recursive = FALSE)
  file_names <- file_list %>% `[[`('name')
  short_descrip %in% file_names
}

drive_folder_exists <- function(target_folder){
  length(which(as_vector(drive_ls()['name']) == target_folder)) != 0
}

ee_to_drive_to_local <- function(ee_object, drive_description, drive_folder, 
                                 drive_fileNamePrefix = NULL, drive_timePrefix = FALSE, 
                                 drive_selectors = NULL,
                                 region = NULL, scale = NULL, maxPixels = NULL, 
                                 local_folder = getwd(), overwrite = FALSE, 
                                 ee_type, # = c('image', 'table'),
                                 file_extension = NULL){
  # ee_type <- match.arg(ee_type)
  dir.create(file.path(local_folder))
  if(nchar(raster::extension(drive_description)) == 0){
    if(is.null(file_extension)){
      stop('file_extension must be specified.\n')
    } else if(file_extension == ''){
      stop('file_extension must be specified.\n')
    } else if (!grepl('[.]', file_extension)){
      file_extension <- paste0(".", file_extension)
    }
    drive_description <- paste0(drive_description, file_extension)
  } else  if(nchar(raster::extension(drive_description)) > 0){
    true_file_extension <- raster::extension(drive_description)
    if (is.null(file_extension)){
      file_extension <- true_file_extension
    } else if (file_extension == ''){
      file_extension <- true_file_extension
    } else if (file_extension != true_file_extension){
      if (!grepl('[.]', file_extension)){
        file_extension <- paste0(".", file_extension)
        warning('File extension and the specified file extension do not match, using the specified file extension.\n')
      }
    }
    drive_description <- paste0(tools::file_path_sans_ext(drive_description), file_extension)
  }
  short_descrip = tail(unlist(strsplit(drive_description, '/')), n =1)
  base_path = gsub(pattern = paste0('/', short_descrip), replacement = '', x = drive_description)
  local_path = paste0(local_folder, '/', short_descrip)
  # return(drive_description)
  if(!drive_file_exists(drive_description)){
    if (ee_type == 'image'){
      task_vector = ee_image_to_drive(
        image = ee_object,
        description = tools::file_path_sans_ext(short_descrip),
        folder = drive_folder,
        fileNamePrefix = drive_fileNamePrefix,
        timePrefix = drive_timePrefix,
        dimensions = NULL,
        region = region,
        scale = scale,
        maxPixels = maxPixels,
        fileFormat = 'GeoTIFF'
      )
      task_vector$start()
      ee_task_monitoring(task_vector)
    } else if (ee_type == 'table'){
      if(raster::extension(drive_description) %in% c('.shp')){
        task_vector = ee_table_to_drive(
          collection = ee_object,
          description = tools::file_path_sans_ext(short_descrip),
          folder = drive_folder,
          fileNamePrefix = drive_fileNamePrefix,
          timePrefix = drive_timePrefix,
          fileFormat = 'SHP',
          selectors = drive_selectors
        )
        task_vector$start()
        ee_task_monitoring(task_vector)
      } else if (raster::extension(drive_description) %in% c('.csv', '.CSV')){
        task_vector = ee_table_to_drive(
          collection = ee_object,
          description = tools::file_path_sans_ext(short_descrip),
          folder = drive_folder,
          fileNamePrefix = drive_fileNamePrefix,
          timePrefix = drive_timePrefix,
          fileFormat = 'CSV',
          selectors = drive_selectors
        )
        task_vector$start()
        ee_task_monitoring(task_vector)
      }
    }
    # task_vector$start()
    # ee_task_monitoring(task_vector)
  }
  # task_vector$start()
  # ee_task_monitoring(task_vector)
  # local_path
  if(!file.exists(local_path)){
    drive_download(
      file = drive_description,
      path = local_path
    )
  } else {
    if(!overwrite){
      warning('file exist, use overwrite = TRUE to overwrite it.\n')
    }
  }
}

gplot_data <- function(x, maxpixels=ncell(x),...){
  if (maxpixels < ncell(x)){
    x <- sampleRegular(x, size = maxpixels, asRaster=TRUE)
  }
  coords <- xyFromCell(x, seq_len(ncell(x)))
  dat <- stack(as.data.frame(getValues(x)))
  names(dat) <- c('value', 'variable')
  
  dat <- cbind(coords, dat)
}

