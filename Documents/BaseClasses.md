# Base Class Method Reference

The three base classes form a hierarchy.  All public methods are inherited by
the concrete product classes unless overridden.

```
nisarBaseHDF
├── nisarBaseRangeDopplerHDF  (adds RD geometry + geolocation)
└── nisarBaseGeocodedHDF      (adds GEO grid + projection)
```

---

## nisarBaseHDF

**Source:** `nisarBaseHDF.py`

The abstract root class.  Handles HDF5 I/O, metadata parsing, data loading,
writing, masking, ionosphere correction, and SAR geometry calculations that
are common to all product types.

### Opening

```python
openHDF(hdfFile, noLoadData=False, fields=None,
        useRos3=False, page_buf_size=2*1024**3,
        downsampleFactor={'downsampleFactorRow':1,'downsampleFactorColumn':1},
        **product_keywords)
```

| Parameter | Description |
|---|---|
| `hdfFile` | Local path, `s3://`, or `https://` URL |
| `noLoadData` | Parse metadata only; do not load data arrays |
| `fields` | List of field names to load (None = product defaults) |
| `useRos3` | Use HDF5 ROS3 driver for S3 (slower but less memory) |
| `page_buf_size` | HDF5 page buffer (bytes); set to file size for best perf |
| `downsampleFactor` | Dict `{downsampleFactorRow, downsampleFactorColumn}` or int |

Remote files are downloaded to memory via `_openS3InMemory` or
`_openHTTPInMemory` unless `useRos3=True`.

### Writing

```python
writeData(filenameRoot, bands=None,
          tiff=True, byteOrder='LSB', grimp=False,
          noSuffix=False, driverName='COG',
          quickLook=False, scale=1.,
          vrtOnly=False, vrtFile=None,
          scaleToPixels=False,   # ROFF/GOFF only
          layers=None)           # ROFF/GOFF only
```

| Parameter | Description |
|---|---|
| `filenameRoot` | Base name; band suffix and extension appended |
| `bands` | Fields to write (None = all loaded fields) |
| `tiff` | Write GeoTIFF if True, flat binary if False |
| `driverName` | `'COG'` (default) or `'GTiff'` |
| `byteOrder` | `'LSB'` (default) or `'MSB'` for binary output |
| `grimp` | Shift GEO origin to lower-left (GrIMP binary convention) |
| `noSuffix` | Omit band name from filename (single-band output only) |
| `quickLook` | Write PNG quick-look instead |
| `scale` | Multiply data by this factor before writing |
| `vrtOnly` | Write VRT pointing at HDF5 bands; no data extracted |
| `vrtFile` | Override default VRT filename (`filenameRoot.vrt`) |
| `scaleToPixels` | Convert offsets to pixels (ROFF/GOFF) |
| `layers` | Layer numbers to write (ROFF/GOFF, default all loaded) |

A `.vrt` is always written alongside the data files (or alone when `vrtOnly=True`).

### Data loading and refreshing

```python
# Reload a subset of fields (e.g. after changing polarization)
obj.refresh(dataFields=['unwrappedPhase', 'coherenceMagnitude'])

# Get image data for a single field into self.fieldName
obj.getImageData('unwrappedPhase', useNumpy=True)
```

### Masking and corrections

```python
maskPhase(largest=True)
# Zero out phase where connectedComponents == 0.
# largest=True: only keep the largest connected component.

applyMask(maskFile)
# Apply an external binary raster mask (any GDAL-readable format).

ionosphereCorrectPhase(masked=False)
# Subtract ionospherePhaseScreen from unwrappedPhase (or wrappedInterferogram).
# masked=True: apply coherence mask first.

removeOutlierOffsets(filterField, threshold, ...)
# Filter offset fields by a quality metric threshold.
```

### Metadata helpers

```python
printParams()                  # Pretty-print all loaded parameters
assembleMeta(bands=None)       # Build self.meta dict for VRT metadata
readGeodatAsDataFrame()        # Read .geodat file as pandas DataFrame
```

### Coordinate conversions

```python
lltoxy(lat, lon)   # lat/lon (°) → projected x/y (km)
xytoll(x, y)       # projected x/y (km) → lat/lon (°)
```

### SAR geometry (orbit-based)

```python
# Compute look/incidence angles at given range, time, height
lookAngle, incAngle = obj.computeAngles(slantRange, zeroDopplerTime, z)

# Satellite position and velocity at a given zero-Doppler time
pos, vel = obj.getSatPositionAndVel(zdTime)

# Scene centre satellite height (km above WGS84)
obj.getSceneCenterSatelliteHeight()

# Range-time → ECEF, then to lat/lon/height
lat, lon, h = obj.RTtoLatLon(slantRange, zdTime, z)
```

---

## nisarBaseRangeDopplerHDF

**Source:** `nisarBaseRangeDopplerHDF.py`  
Inherits from `nisarBaseHDF`. Adds geometry specific to range-Doppler
(slant-range / zero-Doppler time) coordinate products.

### Geotransform

```python
getGeoTransform(tiff=True, grimp=True)
# Returns [rOrigin, dR, 0, aOrigin, 0, dA]
# tiff=False: dA > 0 (azimuth time increases downward — GrIMP binary convention)
# tiff=True:  dA < 0 (north-up GDAL convention, origin at last azimuth time)
```

### Grid setup and interpolation

```python
setupRangeDopplerGrid()
# Create 1-D range and zero-Doppler time vectors for the ML grid.

setupDataCube(cubeName, zGridPts=None)
# Load and set up a 3-D HDF5 metadata cube (e.g. 'incidenceAngle').

interpGrid(interpolator, slantRange, zeroDopplerTime, z)
# Trilinearly interpolate a metadata cube at (range, time, height) points.

xyCube(slantRange, zeroDopplerTime, z, image='reference')
# Map range/time/height → projected x/y via orbit geometry.
```

### Convenience cube methods

All accept arrays of (x, y, z) **or** (slantRange, zeroDopplerTime, z) points
and return interpolated arrays:

| Method | Output |
|---|---|
| `incidenceAngleCube(x, y, z)` | Incidence angle (rad) |
| `elevationAngleCube(x, y, z)` | Elevation angle (rad) |
| `losUnitVectorCube(x, y, z)` | LOS unit vector [3, …] |
| `alongTrackUnitVectorCube(x, y, z)` | Along-track unit vector [3, …] |
| `baselineCube(x, y, z)` | Perpendicular baseline (m) |
| `groundTrackVelocityCube(x, y, z)` | Ground-track velocity (m/s) |

### Geolocation

```python
RDtoLatLon(slantRange, zdTime, z)
# Slant-range + zero-Doppler time + height → lat, lon (degrees)

getCorners()        # Corner lat/lons of the scene
getCenterLatLon()   # Scene centre lat/lon

getGeojson()        # Build GeoJSON feature for scene footprint
writeGeodatGeojson(filename)  # Write GrIMP .geojson geodat file
```

### Timing corrections

```python
getCorrectedTime()       # Apply any zero-Doppler time corrections
getRangeErrorCorrection()  # Compute range timing correction
getSkewOffsets()         # Compute geometric skew offset
```

---

## nisarBaseGeocodedHDF

**Source:** `nisarBaseGeocodedHDF.py`  
Inherits from `nisarBaseHDF`. Adds geometry for geocoded (projected x/y) products.

### Geotransform

```python
getGeoTransform(tiff=True, grimp=False)
# Returns [x0, dx, 0, y0, 0, dy]
# tiff=True (default): dy < 0 (north-up)
# grimp=True: origin at lower-left corner (GrIMP binary convention)
```

### Grid setup and interpolation

```python
setupXYGrid()
# Build 1-D x and y coordinate vectors.

setupDataCube(cubeName, zGridPts=None)
# Load and prepare a 3-D HDF5 metadata cube.

interpGrid(interpolator, xGrid, yGrid, z)
# Trilinearly interpolate at (x, y, height) points.

rangeDopplerCube(x, y, z, image='reference')
# Map projected x/y/height → slant range + zero-Doppler time via orbit.
```

### Convenience cube methods

Same as the range-Doppler variants; accept (x, y, z) in the projected CRS:

`incidenceAngleCube`, `elevationAngleCube`, `losUnitVectorCube`,
`alongTrackUnitVectorCube`, `baselineCube`, `groundTrackVelocityCube`

### Coordinate queries

```python
getExtent()         # (xmin, xmax, ymin, ymax) of the geocoded grid (m)
getGeoCoordinates() # Arrays of x and y coordinates for every pixel
writeGeodat(filename)  # Write GrIMP .geodat binary metadata file
```
