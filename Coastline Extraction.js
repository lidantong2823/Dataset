/************************************************************
 * Multi-Sensor Landsat Coastline Extraction
 *
 * Sensor strategy:
 * 2000–2012: Landsat 7 ETM+ primary + Landsat 5 TM gap filling
 * 2013:      Landsat 7 ETM+ only
 * 2014+:     Landsat 8 OLI
 *
 * Workflow:
 * Landsat preprocessing -> annual median composite ->
 * L5 gap filling -> MNDWI/NDVI/DEM filtering ->
 * morphological filtering -> Canny coastline detection ->
 * raster/vector export
 ************************************************************/


// 1. Configuration
var CONFIG = {
  year: 2015,
  roiAsset: 'users/your_username/your_roi_asset',
  maxCloudCover: 30,
  mndwiThreshold: 0.25,
  ndviThreshold: 0.10,
  maxElevation: 20,
  minLandPatchPixels: 25,
  smoothingRadiusMeters: 30,
  holeFillRadiusPixels: 5,
  cannyThreshold: 0.5,
  cannySigma: 1,
  scale: 30,
  exportCrs: 'EPSG:4326',
  exportFolder: 'GEE_Coastline_Export'
};


// 2. Study Area and Date Range
var roi = ee.FeatureCollection(CONFIG.roiAsset);
var roiGeometry = roi.geometry();
var startDate = ee.Date.fromYMD(CONFIG.year, 1, 1);
var endDate = startDate.advance(1, 'year');

Map.centerObject(roi, 9);
Map.addLayer(roi, {color: 'gray'}, 'Study Area', false);


// 3. Common Band Names
var COMMON_BANDS = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1'];
var ANALYSIS_BANDS = ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'MNDWI', 'NDVI'];


// 4. Landsat 5/7 Preprocessing
function preprocessL57(image) {
  var qa = image.select('QA_PIXEL');
  var qaMask = qa.bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);
  var opticalBands = image.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5'])
    .multiply(0.0000275).add(-0.2).rename(COMMON_BANDS);

  return opticalBands.updateMask(qaMask).updateMask(saturationMask)
    .copyProperties(image, image.propertyNames());
}


// 5. Landsat 8 Preprocessing
function preprocessL8(image) {
  var qa = image.select('QA_PIXEL');
  var qaMask = qa.bitwiseAnd(parseInt('11111', 2)).eq(0);
  var saturationMask = image.select('QA_RADSAT').eq(0);
  var opticalBands = image.select(['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6'])
    .multiply(0.0000275).add(-0.2).rename(COMMON_BANDS);

  return opticalBands.updateMask(qaMask).updateMask(saturationMask)
    .copyProperties(image, image.propertyNames());
}


// 6. Spectral Indices
function addSpectralIndices(image) {
  var mndwi = image.normalizedDifference(['Green', 'SWIR1']).rename('MNDWI');
  var ndvi = image.normalizedDifference(['NIR', 'Red']).rename('NDVI');
  return image.addBands([mndwi, ndvi]);
}


// 7. Empty Fallback Composite
function createEmptyComposite() {
  return ee.Image.constant([0, 0, 0, 0, 0, 0, 0]).rename(ANALYSIS_BANDS)
    .updateMask(ee.Image.constant(0));
}


// 8. Landsat Collections
function getLandsat7Collection() {
  return ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
    .filterBounds(roiGeometry)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', CONFIG.maxCloudCover))
    .map(preprocessL57)
    .map(addSpectralIndices);
}

function getLandsat5Collection() {
  return ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
    .filterBounds(roiGeometry)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', CONFIG.maxCloudCover))
    .map(preprocessL57)
    .map(addSpectralIndices);
}

function getLandsat8Collection() {
  return ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
    .filterBounds(roiGeometry)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt('CLOUD_COVER', CONFIG.maxCloudCover))
    .map(preprocessL8)
    .map(addSpectralIndices);
}


// 9. Annual Composite
function createAnnualComposite(year) {
  if (year <= 2013) {
    var l7Collection = getLandsat7Collection();
    var l7Median = l7Collection.median();
    print('Landsat 7 scene count:', l7Collection.size());

    if (year <= 2012) {
      var l5Collection = getLandsat5Collection();
      var l5Median = ee.Image(ee.Algorithms.If(
        l5Collection.size().gt(0),
        l5Collection.median(),
        createEmptyComposite()
      ));

      print('Landsat 5 gap-fill scene count:', l5Collection.size());

      return l7Median.unmask(l5Median).set({
        primary_sensor: 'Landsat 7 ETM+',
        gap_fill_sensor: 'Landsat 5 TM',
        year: year
      });
    }

    return l7Median.set({
      primary_sensor: 'Landsat 7 ETM+',
      gap_fill_sensor: 'None',
      year: year
    });
  }

  var l8Collection = getLandsat8Collection();
  print('Landsat 8 scene count:', l8Collection.size());

  return l8Collection.median().set({
    primary_sensor: 'Landsat 8 OLI',
    gap_fill_sensor: 'None',
    year: year
  });
}


// 10. Generate Annual Composite
var annualComposite = createAnnualComposite(CONFIG.year);
var medianMndwi = annualComposite.select('MNDWI');
var medianNdvi = annualComposite.select('NDVI');

print('Processing year:', CONFIG.year);
print('Annual composite:', annualComposite);


// 11. DEM and Water Extraction
var dem = ee.Image('USGS/SRTMGL1_003').select('elevation');

var rawWaterMask = medianMndwi.gt(CONFIG.mndwiThreshold)
  .and(medianNdvi.lt(CONFIG.ndviThreshold))
  .and(dem.lt(CONFIG.maxElevation));


// 12. Spatial Noise Filtering
var landMask = rawWaterMask.not();
var landPatchSize = landMask.connectedPixelCount(CONFIG.minLandPatchPixels, true);
var cleanedLandMask = landMask.where(landPatchSize.lt(CONFIG.minLandPatchPixels), 0);
var waterMask = cleanedLandMask.not().focal_mode(CONFIG.smoothingRadiusMeters, 'circle', 'meters');


// 13. Morphological Hole Filling

var filledLandMask = waterMask.not()
  .focal_max({radius: CONFIG.holeFillRadiusPixels, kernelType: 'square', units: 'pixels'})
  .focal_min({radius: CONFIG.holeFillRadiusPixels, kernelType: 'square', units: 'pixels'});

var filledWaterMask = filledLandMask.not();

=
// 14. Canny Coastline Detection
var cannyEdges = ee.Algorithms.CannyEdgeDetector({
  image: filledWaterMask,
  threshold: CONFIG.cannyThreshold,
  sigma: CONFIG.cannySigma
});

var coastlineRaster = cannyEdges.gt(0).selfMask().clip(roiGeometry);


// 15. Vectorization

var waterPolygons = waterMask.selfMask().toByte().reduceToVectors({
  geometry: roiGeometry,
  crs: CONFIG.exportCrs,
  scale: CONFIG.scale,
  geometryType: 'polygon',
  eightConnected: true,
  labelProperty: 'water_class',
  maxPixels: 1e13
});

var islandPolygons = filledLandMask.selfMask().toByte().reduceToVectors({
  geometry: roiGeometry,
  crs: CONFIG.exportCrs,
  scale: CONFIG.scale,
  geometryType: 'polygon',
  eightConnected: false,
  labelProperty: 'land_class',
  maxPixels: 1e13
});


// 16. Visualization
var trueColor = annualComposite.select(['Red', 'Green', 'Blue']);

Map.addLayer(trueColor.clip(roiGeometry), {min: 0.0, max: 0.3, gamma: 1.4}, CONFIG.year + ' Annual Landsat Composite', false);
Map.addLayer(waterMask.selfMask().clip(roiGeometry), {palette: ['0000FF']}, CONFIG.year + ' Water Mask', false);
Map.addLayer(coastlineRaster, {min: 0, max: 1, palette: ['FF0000']}, CONFIG.year + ' Coastline');
Map.addLayer(filledLandMask.selfMask().clip(roiGeometry), {palette: ['00FF00']}, CONFIG.year + ' Filled Land Mask', false);


// 17. Export
var coastlineExport = coastlineRaster.unmask(0).toByte();

Export.image.toDrive({
  image: coastlineExport,
  description: 'Coastline_Raster_' + CONFIG.year,
  folder: CONFIG.exportFolder,
  region: roiGeometry,
  scale: CONFIG.scale,
  maxPixels: 1e13,
  crs: CONFIG.exportCrs
});

Export.table.toDrive({
  collection: waterPolygons,
  description: 'Water_Polygons_' + CONFIG.year,
  folder: CONFIG.exportFolder,
  fileFormat: 'SHP'
});

Export.table.toDrive({
  collection: islandPolygons,
  description: 'Island_Polygons_' + CONFIG.year,
  folder: CONFIG.exportFolder,
  fileFormat: 'SHP'
});