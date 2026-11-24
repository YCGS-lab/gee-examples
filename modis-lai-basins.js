// Point near the middle of CT
var ct_center = ee.Geometry.Point([-72.66865061353076, 41.76172521301823]);

// New England-ish bounding box
var b_box = ee.Geometry.Polygon(
        [[[-75.46892570546962, 44.24404942901968],
          [-75.46892570546962, 40.3369472454705],
          [-69.39348625234462, 40.3369472454705],
          [-69.39348625234462, 44.24404942901968]]], null, false);

// Load USGS basins as an example
var huc10 = ee.FeatureCollection('USGS/WBD/2017/HUC10');

// MODIS 4-day Leaf Area Index (LAI)
var modis_lai = ee.ImageCollection("MODIS/061/MCD15A3H");

// Subset to New England-ish region
var huc10_sub = huc10.filterBounds(b_box);

// Add a sequential integer index (1, 2, 3, ...) to each HUC code (e.g., 0106000104).
// This lets us construct a safe and compact image in the `paint` step,
// without worrying about the full HUC code getting mangled.
// We will retrieve the original HUC code at the end.
var hucCodes = huc10_sub.aggregate_array('huc10');
var basins = huc10_sub.map(function(f) {
  return f.set('zone', hucCodes.indexOf(f.get('huc10')));
});

// Convert HUC10 zones to a raster image
var zones = ee.Image().int16().paint(basins, 'zone').rename('zone');

// Get MODIS's native projection
var modisProj = modis_lai.first().select('Lai').projection();

// Get the complete encompassing region of the HUC10 subset
// (note: not b_box, because we include basins that partially intersect our bounding box)
var region = huc10_sub.geometry().bounds();

// Construct the "reducers", which will aggregate the results by basin.
// We want three values here: lai_mean, fpar_mean, and n_pixels (because lai and fpar are masked together).
// First, we construct the reducers separately.
var lai_mean = ee.Reducer.mean().setOutputs(['lai_mean']);
var counts = ee.Reducer.count().setOutputs(['n_pixels']);
var fpar_mean = ee.Reducer.mean().setOutputs(['fpar_mean']);

// Now, combine them into a single big reducer.
var reducer = lai_mean
// NOTE: We need to combine lai_mean and counts first because they have a shared input (hence: sharedInputs: true)
  .combine({reducer2: counts, sharedInputs: true})
// After this step, the reducer takes one shared input and returns 2 values.
// NOTE: We add fpar next because it uses a separate input
  .combine({reducer2: fpar_mean})
// After this step, the reducer takes *two* different inputs; one shared input for LAI and counts, and a second for fPAR.
// The final step applies the reducer by group. Here, we specify that the grouping variable will be at index 2 (3rd item).
  .group({groupField: 2, groupName: 'zone'});
  
print(reducer.getOutputs());
  
// Function that takes a single MODIS image (one timestep) and returns a
// featureCollection ("shapefile") of the mean LAI, fPAR, and number of pixels
// for every HUC basin.
var calcBasinStats = function (img) {
  var qc = img.select('FparLai_QC');
  var lai = img.select('Lai');
  var fpar = img.select('Fpar');
  
  // Screen only valid values
  var valid = qc.bitwiseAnd(1).eq(0)
  // Values greater than 100 for LAI and fPAR are fill codes, not real data.
    .and(lai.lte(100))
    .and(fpar.lte(100));
    
  var stack = lai.multiply(0.1).rename('lai')
    .addBands(fpar.multiply(0.01).rename('fpar'))
    .addBands(zones)
    // NOTE: Apply the mask at the end.
    // Otherwise, you can end up with blank stacks that error the `get` step later.
    .updateMask(valid);
    
  var stats = stack.reduceRegion({
    reducer: reducer,
    geometry: region,
    crs: modisProj,
    maxPixels: 1e10
  });
  
  var date = img.date().format('YYYY-MM-dd');
  
  return ee.FeatureCollection(ee.List(stats.get('groups')).map(function (g) {
    g = ee.Dictionary(g);
    var lai = g.get('lai_mean', null);
    var fpar = g.get('fpar_mean', null);
    var npix = g.get('n_pixels', 0);
    return ee.Feature(null, {
      huc10: hucCodes.get(ee.Number(g.get('zone')).int()),
      date: date,
      lai: lai,
      fpar: fpar,
      n_valid: npix
    });
  }));
};

// Quick sanity check on a short window before exporting everything
// print(ee.FeatureCollection(
//   modis_lai.filterDate('2020-07-01', '2020-07-10').map(calcBasinStats)
// ).flatten().limit(10));
  
var ts = ee.FeatureCollection(modis_lai.filterDate('2020-01-01', '2020-07-01').map(calcBasinStats)).flatten();

Export.table.toDrive({
  collection: ts,
  description: 'huc10_modis_lai_fpar',
  fileFormat: 'CSV',
  selectors: ['huc10', 'date', 'lai', 'fpar', 'n_valid']
});

// 429 basins x 6 months of MODIS data = ~38 sec runtime /1153 EECU-seconds.
