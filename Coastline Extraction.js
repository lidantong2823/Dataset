/************************************************************
 * Multi-Sensor Landsat Island Boundary and Coastline Extraction
 * Google Earth Engine JavaScript API
 *
 * Sensor strategy:
 * 2000-2002: Landsat 7 ETM+
 * 2003-2012: Landsat 7 ETM+ with Landsat 5 TM gap filling
 * 2013:      Landsat 7 ETM+
 * 2014-2025: Landsat 8 OLI
 *
 * Workflow:
 * Annual reflectance composite -> MNDWI + AWEI ->
 * Canny-guided Otsu thresholds -> combined water mask ->
 * morphological closing -> coastline raster + island polygons
 ************************************************************/

// ==========================================================
// 1. Configuration
// ==========================================================
var CONFIG = {
  year: 2020,
  roiAsset: 'users/your_username/your_roi_asset',
  usePercentileComposite: true,
  compositePercentile: 10,
  gapFillKernelPixels: 10,
  gapFillWindowDays: 10,
  gapFillMinNeighbors: 144,
  gapFillMinScale: 1 / 3,
  gapFillMaxScale: 3,
  otsuCannyThreshold: 0.7,
  otsuCannySigma: 1,
  coastlineCannyThreshold: 0.5,
  coastlineCannySigma: 1,
  closingRadiusMeters: 30,
  scale: 30,
  exportCrs: 'EPSG:4326',
  exportFolder: 'GEE_Coastline_Export'
};

if (CONFIG.year < 2000 || CONFIG.year > 2025) {
  throw new Error('CONFIG.year must be between 2000 and 2025.');
}

// ==========================================================
// 2. Study Area
// ==========================================================
var roi = ee.FeatureCollection(CONFIG.roiAsset);
var region = roi.geometry();

Map.centerObject(roi, 8);
Map.addLayer(roi, {color: 'FFFF00'}, 'Study Area', false);

// ==========================================================
// 3. Landsat Preprocessing
// ==========================================================
function preprocessL57(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);
  var optical = image.select(['SR_B2', 'SR_B4', 'SR_B5', 'SR_B7'])
    .multiply(0.0000275).add(-0.2)
    .rename(['green', 'nir', 'swir1', 'swir2']);

  return optical.updateMask(qaMask).updateMask(saturationMask)
    .copyProperties(image, image.propertyNames());
}

function preprocessL8(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);
  var optical = image.select(['SR_B3', 'SR_B5', 'SR_B6', 'SR_B7'])
    .multiply(0.0000275).add(-0.2)
    .rename(['green', 'nir', 'swir1', 'swir2']);

  return optical.updateMask(qaMask).updateMask(saturationMask)
    .copyProperties(image, image.propertyNames());
}

// ==========================================================
// 4. Landsat 7 Gap Filling with Landsat 5
// ==========================================================
function gapFillL7WithL5(sourceCollection, fillCollection) {
  var bands = ['green', 'nir', 'swir1', 'swir2'];
  var kernel = ee.Kernel.square(CONFIG.gapFillKernelPixels * CONFIG.scale, 'meters', false);

  return sourceCollection.map(function(image) {
    var date = ee.Date(image.get('system:time_start'));

    var sourceWindow = sourceCollection
      .filterBounds(region)
      .filterDate(date.advance(-CONFIG.gapFillWindowDays, 'day'), date.advance(CONFIG.gapFillWindowDays, 'day'));

    var fillWindow = fillCollection
      .filterBounds(region)
      .filterDate(date.advance(-CONFIG.gapFillWindowDays, 'day'), date.advance(CONFIG.gapFillWindowDays, 'day'));

    var canFill = sourceWindow.size().gt(0).and(fillWindow.size().gt(0));

    var filled = ee.Image(ee.Algorithms.If(canFill, (function() {
      var sourceImage = sourceWindow.median().select(bands);
      var fillImage = fillWindow.median().select(bands);
      var commonMask = sourceImage.mask().and(fillImage.mask());

      var sourceMasked = sourceImage.updateMask(commonMask);
      var fillMasked = fillImage.updateMask(commonMask);

      var regression = fillMasked.addBands(sourceMasked);
      regression = regression.select(regression.bandNames().sort());

      var fit = regression.reduceNeighborhood(
        ee.Reducer.linearFit().forEach(sourceImage.bandNames()),
        kernel,
        null,
        false
      );

      var offset = fit.select('.*_offset');
      var scale = fit.select('.*_scale');

      var statsReducer = ee.Reducer.mean().combine(ee.Reducer.stdDev(), null, true);
      var sourceStats = sourceImage.reduceNeighborhood(statsReducer, kernel, null, false);
      var fillStats = fillImage.reduceNeighborhood(statsReducer, kernel, null, false);

      var scaleFactor = sourceStats.select('.*stdDev').divide(fillStats.select('.*stdDev'));
      var offsetFactor = sourceStats.select('.*mean')
        .subtract(fillStats.select('.*mean').multiply(scaleFactor));

      var invalidScale = scale.lt(CONFIG.gapFillMinScale).or(scale.gt(CONFIG.gapFillMaxScale));
      scale = scale.where(invalidScale, scaleFactor);
      offset = offset.where(invalidScale, offsetFactor);

      var fallbackScale = scale.lt(CONFIG.gapFillMinScale).or(scale.gt(CONFIG.gapFillMaxScale));
      scale = scale.where(fallbackScale, 1);
      offset = offset.where(
        fallbackScale,
        sourceStats.select('.*mean').subtract(fillStats.select('.*mean'))
      );

      var neighborCount = commonMask.reduceNeighborhood(
        ee.Reducer.count(),
        kernel,
        null,
        true,
        'boxcar'
      );

      var scaledFill = fillImage.multiply(scale).add(offset)
        .updateMask(neighborCount.gte(CONFIG.gapFillMinNeighbors));

      return image.unmask(scaledFill, true);
    })(), image));

    return filled.copyProperties(image, ['system:time_start']);
  });
}

// ==========================================================
// 5. Annual Sensor Selection
// ==========================================================
function getAnnualCollection(year) {
  var start = year + '-01-01';
  var end = (year + 1) + '-01-01';

  if (year >= 2000 && year <= 2002) {
    return {
      sensor: 'Landsat7_ETM+',
      collection: ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
        .filterBounds(region)
        .filterDate(start, end)
        .map(preprocessL57)
    };
  }

  if (year >= 2003 && year <= 2012) {
    var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
      .filterBounds(region)
      .filterDate(start, end)
      .map(preprocessL57);

    var l5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
      .filterBounds(region)
      .filterDate(start, end)
      .map(preprocessL57);

    return {
      sensor: 'Landsat7_ETM+_L5GapFill',
      collection: gapFillL7WithL5(l7, l5)
    };
  }

  if (year === 2013) {
    return {
      sensor: 'Landsat7_ETM+',
      collection: ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
        .filterBounds(region)
        .filterDate(start, end)
        .map(preprocessL57)
    };
  }

  return {
    sensor: 'Landsat8_OLI',
    collection: ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
      .filterBounds(region)
      .filterDate(start, end)
      .map(preprocessL8)
  };
}

// ==========================================================
// 6. Annual Composite
// ==========================================================
function makeAnnualComposite(collection) {
  var bands = ['green', 'nir', 'swir1', 'swir2'];

  if (CONFIG.usePercentileComposite) {
    return collection.select(bands)
      .reduce(ee.Reducer.percentile([CONFIG.compositePercentile]))
      .rename(bands)
      .set('composite', 'p' + CONFIG.compositePercentile + '_reflectance');
  }

  return collection.select(bands).median().rename(bands)
    .set('composite', 'annual_median_reflectance');
}

// ==========================================================
// 7. Water Indices
// ==========================================================
function addWaterIndices(image) {
  var mndwi = image.normalizedDifference(['green', 'swir1']).rename('mndwi');

  var awei = image.expression(
    '4 * (GREEN - SWIR1) - (0.25 * NIR + 2.75 * SWIR2)',
    {
      GREEN: image.select('green'),
      NIR: image.select('nir'),
      SWIR1: image.select('swir1'),
      SWIR2: image.select('swir2')
    }
  ).rename('awei');

  return image.addBands([mndwi, awei]);
}

// ==========================================================
// 8. Otsu Thresholding
// ==========================================================
function otsu(histogram) {
  histogram = ee.Dictionary(histogram);

  var counts = ee.Array(histogram.get('histogram'));
  var means = ee.Array(histogram.get('bucketMeans'));
  var size = means.length().get([0]);
  var indices = ee.List.sequence(1, ee.Number(size).subtract(1));

  var bss = indices.map(function(i) {
    i = ee.Number(i);

    var counts1 = counts.slice(0, 0, i);
    var means1 = means.slice(0, 0, i);
    var counts2 = counts.slice(0, i);
    var means2 = means.slice(0, i);

    var weight1 = ee.Number(counts1.reduce(ee.Reducer.sum(), [0]).get([0]));
    var weight2 = ee.Number(counts2.reduce(ee.Reducer.sum(), [0]).get([0]));

    var mean1 = ee.Number(
      means1.multiply(counts1).reduce(ee.Reducer.sum(), [0]).get([0])
    ).divide(weight1);

    var mean2 = ee.Number(
      means2.multiply(counts2).reduce(ee.Reducer.sum(), [0]).get([0])
    ).divide(weight2);

    return weight1.multiply(weight2).multiply(mean1.subtract(mean2).pow(2));
  });

  var maxBss = ee.List(bss).reduce(ee.Reducer.max());
  var index = ee.List(bss).indexOf(maxBss);

  return ee.Number(means.toList().get(ee.Number(index).add(1)));
}

function getOtsuThreshold(indexImage, bandName) {
  var edges = ee.Algorithms.CannyEdgeDetector(
    indexImage,
    CONFIG.otsuCannyThreshold,
    CONFIG.otsuCannySigma
  ).focalMax(1);

  var histogram = ee.Dictionary(
    indexImage.updateMask(edges).reduceRegion({
      reducer: ee.Reducer.histogram({maxBuckets: 512}),
      geometry: region,
      scale: CONFIG.scale,
      maxPixels: 1e13,
      tileScale: 4
    }).get(bandName)
  );

  return otsu(histogram);
}

function getCombinedWaterMask(mndwiImage, aweiImage) {
  var mndwiThreshold = getOtsuThreshold(mndwiImage, 'mndwi');
  var aweiThreshold = getOtsuThreshold(aweiImage, 'awei');

  var waterMask = mndwiImage.gte(mndwiThreshold)
    .and(aweiImage.gte(aweiThreshold))
    .rename('water')
    .toByte();

  return {
    mndwiThreshold: mndwiThreshold,
    aweiThreshold: aweiThreshold,
    waterMask: waterMask
  };
}

// ==========================================================
// 9. Morphological Closing
// ==========================================================
function morphologicalClosing(mask) {
  return mask
    .focalMax({radius: CONFIG.closingRadiusMeters, units: 'meters'})
    .focalMin({radius: CONFIG.closingRadiusMeters, units: 'meters'})
    .rename('water_closed')
    .toByte();
}

// ==========================================================
// 10. Remove Interior Polygon Holes
// ==========================================================
function removeHoles(geometry) {
  var type = geometry.type();

  return ee.Geometry(ee.Algorithms.If(
    ee.String(type).compareTo('Polygon').eq(0),

    ee.Geometry.Polygon(
      [ee.List(ee.List(geometry.coordinates()).get(0))],
      null,
      false
    ),

    ee.Algorithms.If(
      ee.String(type).compareTo('MultiPolygon').eq(0),

      ee.Geometry.MultiPolygon(
        ee.List(geometry.coordinates()).map(function(polygon) {
          polygon = ee.List(polygon);
          return ee.List([ee.List(polygon.get(0))]);
        }),
        null,
        false
      ),

      geometry
    )
  ));
}

// ==========================================================
// 11. Island Polygon Vectorization
// ==========================================================
function vectorizeLand(mask, year, sensor) {
  var polygons = mask.selfMask().reduceToVectors({
    geometry: region,
    scale: CONFIG.scale,
    geometryType: 'polygon',
    eightConnected: true,
    labelProperty: 'class',
    reducer: ee.Reducer.countEvery(),
    maxPixels: 1e13,
    tileScale: 4
  });

  return polygons.map(function(feature) {
    return ee.Feature(removeHoles(feature.geometry()), {
      year: year,
      sensor: sensor,
      class: 'island',
      comp: CONFIG.usePercentileComposite ? 'p' + CONFIG.compositePercentile : 'median',
      method: 'MNDWI_AWEI_OTSU',
      shoreline: CONFIG.usePercentileComposite ? 'high_water_proxy' : 'median_waterline_proxy'
    });
  });
}

// ==========================================================
// 12. Run Selected Year
// ==========================================================
function runYear(year) {
  var data = getAnnualCollection(year);
  var reflectanceComposite = makeAnnualComposite(data.collection);
  var composite = addWaterIndices(reflectanceComposite);

  var mndwi = composite.select('mndwi');
  var awei = composite.select('awei');

  var result = getCombinedWaterMask(mndwi, awei);
  var waterClosed = morphologicalClosing(result.waterMask);
  var landMask = waterClosed.eq(0).rename('land').toByte();

  var coastline = ee.Algorithms.CannyEdgeDetector({
    image: waterClosed,
    threshold: CONFIG.coastlineCannyThreshold,
    sigma: CONFIG.coastlineCannySigma
  }).gt(0).selfMask().rename('coastline');

  var islandPolygons = vectorizeLand(landMask, year, data.sensor);

  print('Year:', year);
  print('Sensor:', data.sensor);
  print('Available scenes:', data.collection.size());
  print('MNDWI Otsu threshold:', result.mndwiThreshold);
  print('AWEI Otsu threshold:', result.aweiThreshold);

  Map.layers().reset();
  Map.addLayer(roi, {color: 'FFFF00'}, 'Study Area', false);
  Map.addLayer(mndwi, {min: -1, max: 1, palette: ['8C510A', 'F6E8C3', '01665E']}, year + ' MNDWI', false);
  Map.addLayer(awei, {min: -3, max: 3, palette: ['7F0000', 'FCBBA1', '08306B']}, year + ' AWEI', false);
  Map.addLayer(waterClosed.selfMask(), {palette: ['00BFFF']}, year + ' Water Mask', false);
  Map.addLayer(coastline, {palette: ['FF0000']}, year + ' Coastline', true);
  Map.addLayer(
    islandPolygons.style({color: '00FF00', fillColor: '00000000', width: 2}),
    {},
    year + ' Island Polygons',
    true
  );

  Export.image.toDrive({
    image: coastline.unmask(0).toByte(),
    description: 'coastline_raster_' + year,
    folder: CONFIG.exportFolder,
    fileNamePrefix: 'coastline_raster_' + year,
    region: region,
    scale: CONFIG.scale,
    maxPixels: 1e13,
    crs: CONFIG.exportCrs,
    fileFormat: 'GeoTIFF'
  });

  Export.table.toDrive({
    collection: islandPolygons,
    description: 'island_polygons_' + year,
    folder: CONFIG.exportFolder,
    fileNamePrefix: 'island_polygons_' + year,
    fileFormat: 'SHP'
  });
}

// ==========================================================
// 13. Execute
// ==========================================================
runYear(CONFIG.year);