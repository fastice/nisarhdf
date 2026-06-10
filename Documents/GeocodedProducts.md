# Geocoded Product Classes

All three classes inherit from `nisarBaseGeocodedHDF` and share `coordType = 'GEO'`.
Coordinates are **projected x/y** (polar stereographic for Greenland/Antarctica or
UTM elsewhere). The EPSG code is read from the HDF5 and stored as `self.epsg`.

---

## nisarGCOVHDF — Geocoded Covariance

**Source:** `nisarGCOVHDF.py`  
**HDF5 product key:** `GCOV`  
**bands key:** `grids`  
**lookType:** `None` (uses `xSize`/`ySize` rather than `{lookType}RangeSize`)

### Instantiation

```python
myGCOV = nisarhdf.nisarGCOVHDF(frequency='frequencyA')
```

### Opening

```python
myGCOV.openHDF('NISAR_L2_PR_GCOV_…h5',
               fields=['HHHH', 'VVVV'],   # covariance terms to load
               dB=False,                   # convert to dB on load
               sigma0=False)               # sigma0 instead of gamma0
```

### Available fields

Default: `HHHH`, `VVVV`, `HVHV`, `VHVH`  
Additional: `mask`, `numberOfLooks`, `rtcGammaToSigmaFactor`  
(Only terms matching available polarizations will be present.)

Check `myGCOV.covTerms` for the product-specific list.

### Key attributes

| Attribute | Description |
|---|---|
| `HHHH`, `VVVV`, … | float32 arrays — RTC backscatter power |
| `xSize`, `ySize` | Grid dimensions (pixels) |
| `dx`, `dy` | Pixel spacing (m) |
| `epsg` | EPSG code of the projection |
| `NumberRangeLooks` | Range looks used in processing |
| `NumberAzimuthLooks` | Azimuth looks used in processing |
| `covTerms` | List of covariance terms present in HDF5 |

### Notes

- Pass `dB=True` to `openHDF` or `writeData` to output 10·log₁₀(power).
- Pass `sigma0=True` to convert gamma0 → sigma0 using the `rtcGammaToSigmaFactor`.
- GCOV files are often very large; `--downsampleFactor` is useful for quick looks.
- `fields` accepts both polarization labels and explicit covariance-term names.
  Use `--info` / `printParams()` to see what is available.

---

## nisarGUNWHDF — Geocoded Unwrapped Interferogram

**Source:** `nisarGUNWHDF.py`  
**HDF5 product key:** `GUNW`  
**bands key:** `grids`  
**productType:** `unwrappedInterferogram` (default) or `wrappedInterferogram`  
**lookType:** `ML`

### Instantiation

```python
myGUNW = nisarhdf.nisarGUNWHDF(frequency='frequencyA',
                                productType='unwrappedInterferogram')
# or for wrapped:
myGUNW = nisarhdf.nisarGUNWHDF(productType='wrappedInterferogram')
```

### Available fields

**unwrappedInterferogram:** `unwrappedPhase`, `coherenceMagnitude`,
`connectedComponents`, `ionospherePhaseScreen`,
`ionospherePhaseScreenUncertainty`

Default: `unwrappedPhase`

**wrappedInterferogram:** `wrappedInterferogram`, `coherenceMagnitude`

Default: `wrappedInterferogram`

### Key attributes

| Attribute | Description |
|---|---|
| `unwrappedPhase` | float32 geocoded unwrapped phase (rad) |
| `coherenceMagnitude` | float32 coherence [0–1] |
| `connectedComponents` | uint16 connected-component labels |
| `ionospherePhaseScreen` | float32 ionospheric correction (rad) |
| `xSize`, `ySize` | Geocoded grid dimensions |
| `dx`, `dy` | Pixel spacing (m) |
| `MLRangeSize`, `MLAzimuthSize` | Pre-geocoding ML dimensions |

### Processing methods (inherited)

```python
myGUNW.maskPhase()                # zero out pixels with connectedComponents == 0
myGUNW.ionosphereCorrectPhase()   # subtract ionospherePhaseScreen
myGUNW.applyMask('mask.tif')      # apply external raster mask
```

---

## nisarGOFFHDF — Geocoded Pixel Offsets

**Source:** `nisarGOFFHDF.py`  
**HDF5 product key:** `GOFF`  
**bands key:** `grids` / `productType`: `pixelOffsets`  
**lookType:** `Offset`

### Instantiation

```python
myGOFF = nisarhdf.nisarGOFFHDF(frequency='frequencyA',
                                layer='layer1')
```

### Available fields

Default: `slantRangeOffset`, `alongTrackOffset`  
All: `slantRangeOffset`, `slantRangeOffsetVariance`, `alongTrackOffset`,
`alongTrackOffsetVariance`, `crossOffsetVariance`, `correlationSurfacePeak`, `snr`

### Key attributes

| Attribute | Description |
|---|---|
| `xSize`, `ySize` | Grid dimensions |
| `dx`, `dy` | Pixel spacing (m) |
| `epsg` | Projection EPSG |

### Writing offsets

Same suffix convention as ROFF; per-layer VRTs are created.
Pass `scaleToPixels=True` to convert metric offsets to pixels.

---

## Common Geocoded attributes (all three classes)

| Attribute | Description |
|---|---|
| `xSize`, `ySize` | Grid dimensions (pixels) |
| `dx`, `dy` | Pixel spacing (m) |
| `epsg` | EPSG code of output projection |
| `coordType` | Always `'GEO'` |
| `referenceOrbit` | Reference orbit number |
| `track`, `frame` | Track / frame numbers |
| `datetime` | Scene centre `datetime` |
| `centerLat`, `centerLon` | Scene centre lat/lon |
| `LookDirection` | `'left'` or `'right'` |
| `PassType` | `'ascending'` or `'descending'` |
| `Wavelength` | Radar wavelength (m) |

---

## Geotransform convention

`getGeoTransform(tiff=True)` returns `[x0, dx, 0, y0, 0, dy]` where:
- `x0` = western edge x coordinate (m)
- `dx` = pixel width (positive)
- `y0` = northern edge y coordinate (m)
- `dy` = pixel height (negative for north-up convention)

The `grimp=True` variant shifts the origin to the lower-left corner (GrIMP binary
convention) but `dy` sign is unchanged.

---

## Interpolating metadata cubes onto x/y grids

Geocoded products carry HDF5 metadata *cubes* (height-stratified grids of
incidence angle, baseline, LOS/along-track unit vectors, etc.) that can be
interpolated to any set of x/y/z points:

```python
# Suppose you have an array of (x, y) positions in the projected CRS
# and corresponding DEM heights z
inc = myGUNW.incidenceAngleCube(x, y, z)
los = myGUNW.losUnitVectorCube(x, y, z)
base = myGUNW.baselineCube(x, y, z)
```

These methods return arrays of the same shape as the input coordinate arrays.

For reverse-mapping from range-Doppler to geocoded coordinates use
`rangeDopplerCube(x, y, z)` (inherited from `nisarBaseGeocodedHDF`).
