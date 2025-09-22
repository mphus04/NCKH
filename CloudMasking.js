// ==========================
// Hàm mask theo lớp SCL của Sentinel-2 SR
// ==========================
function maskSCL(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(3)   // shadow
    .and(scl.neq(8))      // medium cloud
    .and(scl.neq(9))      // high cloud
    .and(scl.neq(10))     // cirrus
    .and(scl.neq(11));    // snow
  return image.updateMask(mask);
}

// ==========================
// Hàm mask theo Sentinel-2 Cloud Probability
// ==========================
function maskS2Cloudless(image) {
  var cloudProb = image.select('probability');
  var isCloud = cloudProb.gt(50); // nếu xác suất mây > 50%
  return image.updateMask(isCloud.not());
}


// ==========================
// HÀM LẤY ẢNH SENTINEL THEO NĂM (có mask mây)
// ==========================
function getSentinelImage(year, geometry) {
  var numericYear = ee.Number.parse(year);
  var startDate = ee.Date.fromYMD(numericYear, 1, 1);
  var endDate = ee.Date.fromYMD(numericYear, 12, 31);

  // Nếu năm <= 2018 → dùng S2_HARMONIZED với QA60
  if (numericYear.lte(2018)) {
    var s2 = ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
      .filterBounds(geometry)
      .filterDate(startDate, endDate)
      .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20));

    function maskS2clouds(image) {
      var qa = image.select('QA60');
      var cloudBitMask = 1 << 10;
      var cirrusBitMask = 1 << 11;
      var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
        .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
      return image.updateMask(mask)
        .multiply(0.0001)
        .copyProperties(image, ['system:time_start']);
    }

    var masked = s2.map(maskS2clouds);

    return masked.median()
      .select(['B4', 'B3', 'B2'])
      .clip(geometry);
  }

  // Nếu năm > 2018 → dùng S2_SR với Cloud Probability + SCL
  var s2SR = ee.ImageCollection("COPERNICUS/S2_SR")
    .filterBounds(geometry)
    .filterDate(startDate, endDate)
    .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 20))
    .select([
      'B1','B2','B3','B4','B5','B6','B7','B8',
      'B8A','B9','B11','B12','SCL']);

  var s2Clouds = ee.ImageCollection("COPERNICUS/S2_CLOUD_PROBABILITY")
    .filterBounds(geometry)
    .filterDate(startDate, endDate);

  var masked = s2SR.map(function(img) {
    var cloudProbImg = s2Clouds
      .filter(ee.Filter.eq('system:index', img.get('system:index')))
      .first();

    var cloudMask = ee.Algorithms.If(
      cloudProbImg,
      ee.Image(cloudProbImg).select('probability').lt(50),
      ee.Image(1) // giữ lại nếu không có cloud prob
    );

    var sclMask = img.select('SCL')
      .neq(3)
      .and(img.select('SCL').neq(8))
      .and(img.select('SCL').neq(9))
      .and(img.select('SCL').neq(10))
      .and(img.select('SCL').neq(11));

    var combinedMask = ee.Image(cloudMask).and(sclMask);
    return img.updateMask(combinedMask);
  });

  return ee.ImageCollection(masked)
    .median()
    .select(['B4', 'B3', 'B2'])
    .multiply(0.0001)
    .clip(geometry);
}


//export hàm ra để sử dụng như module
exports.getSentinelImage = getSentinelImage;
