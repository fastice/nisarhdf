# Range-Doppler Product Classes

All four classes inherit from `nisarBaseRangeDopplerHDF` and share `coordType = 'RD'`.
Coordinates are **slant range (m)** and **zero-Doppler time (s of day)**.

---

## nisarRSLCHDF — Single Look Complex

**Source:** `nisarRSLCHDF.py`  
**HDF5 product key:** `RSLC`  
**bands key:** `swaths`  
**lookType:** `SLC`  

### Instantiation

```python
myRSLC = nisarhdf.nisarRSLCHDF(frequency='frequencyA')   # default
```

### Opening

```python
myRSLC.openHDF('NISAR_L1_PR_RSLC_…h5',
               fields=['HH'],           # subset; default = all available pols
               noLoadData=False)        # True = metadata only, no arrays loaded
```

### Available fields / polarizations

`HH`, `VV`, `HV`, `VH` (product-dependent; check `myRSLC.polarizations`).

### Key attributes after `openHDF`

| Attribute | Description |
|---|---|
| `HH` (or `VV`, etc.) | `complex64` numpy array `[azimuth, range]` |
| `SLCRangeSize` | Number of range pixels |
| `SLCAzimuthSize` | Number of azimuth lines |
| `SLCRangePixelSize` | Range pixel spacing (m) |
| `SLCAzimuthPixelSize` | Azimuth pixel spacing (m along-track) |
| `SLCNearRange` | Near-range slant range (m) |
| `SLCFarRange` | Far-range slant range (m) |
| `SLCFirstZeroDopplerTime` | First zero-Doppler time (s of day) |
| `SLCLastZeroDopplerTime` | Last zero-Doppler time (s of day) |
| `SLCZeroDopplerTimeSpacing` | Zero-Doppler time step (s) = 1/PRF |
| `PRF` | Pulse repetition frequency (Hz) |
| `NumberRangeLooks` | Always 1 for SLC |
| `NumberAzimuthLooks` | Always 1 for SLC |
| `polarizations` | List of available pols, e.g. `['HH', 'VV']` |
| `Wavelength` | Radar wavelength (m) |

### Notes

- Data type is `complex64` (CFloat32). Modules that need power can call
  `openHDF(..., power=True)` to compute `abs²` on load.
- `writeSLC(filename)` writes a flat binary `>c8` file (legacy GrIMP method).
- `writeOldParFile(filename)` writes a GAMMA-style parameter file.

---

## nisarRUNWHDF — Unwrapped Interferogram

**Source:** `nisarRUNWHDF.py`  
**HDF5 product key:** `RUNW`  
**bands key:** `swaths` / `productType`: `interferogram`  
**lookType:** `ML`

### Instantiation

```python
myRUNW = nisarhdf.nisarRUNWHDF(frequency='frequencyA',
                                polarization=None)   # None → first like-pol
```

### Available fields

Default: `unwrappedPhase`  
All: `unwrappedPhase`, `coherenceMagnitude`, `connectedComponents`,
`ionospherePhaseScreen`, `ionospherePhaseScreenUncertainty`, `digitalElevationModel`

### Key attributes

| Attribute | Type | Description |
|---|---|---|
| `unwrappedPhase` | float32 | Unwrapped interferometric phase (rad) |
| `coherenceMagnitude` | float32 | Coherence [0–1] |
| `connectedComponents` | uint16 | Connected-component label map |
| `ionospherePhaseScreen` | float32 | Ionospheric phase correction (rad) |
| `digitalElevationModel` | float32 | DEM used in processing (m) |
| `MLRangeSize` | int | Multi-looked range dimension |
| `MLAzimuthSize` | int | Multi-looked azimuth dimension |
| `NumberRangeLooks` | int | Range looks |
| `NumberAzimuthLooks` | int | Azimuth looks |

### Processing methods

```python
myRUNW.maskPhase()                # zero out pixels with connectedComponents == 0
myRUNW.ionosphereCorrectPhase()   # subtract ionospherePhaseScreen
myRUNW.applyMask('mask.tif')      # apply external raster mask
myRUNW.removeOutlierOffsets(...)  # filter phase by coherence threshold
```

### Metadata cube interpolation

```python
# Incidence angle interpolated to the ML data grid
myRUNW.incidenceAngleCube(x, y, z)   # z = DEM values at (x,y) points

# LOS unit vector, along-track unit vector, baseline, etc.
myRUNW.losUnitVectorCube(x, y, z)
myRUNW.baselineCube(x, y, z)
myRUNW.groundTrackVelocityCube(x, y, z)
```

---

## nisarRIFGHDF — Range Interferogram (wrapped)

**Source:** `nisarRIFGHDF.py`  
**HDF5 product key:** `RIFG`  
**bands key:** `swaths` / `productType`: `interferogram`  
**lookType:** `ML`

Identical interface to RUNW; available fields are `wrappedInterferogram` and
`coherenceMagnitude`. Shares the same base-class processing methods.

### Available fields

Default: `wrappedInterferogram`  
All: `wrappedInterferogram`, `coherenceMagnitude`

---

## nisarROFFHDF — Range-Doppler Pixel Offsets

**Source:** `nisarROFFHDF.py`  
**HDF5 product key:** `ROFF`  
**bands key:** `swaths` / `productType`: `pixelOffsets`  
**lookType:** `Offset`

### Instantiation

```python
myROFF = nisarhdf.nisarROFFHDF(frequency='frequencyA',
                                layer='layer1')   # layer1–3 available
```

`openHDF` accepts `layers=['layer1','layer2','layer3']` to load multiple layers
at once; each field becomes a 3-D array `[nLayers, azimuth, range]`.

### Available fields

Default: `slantRangeOffset`, `alongTrackOffset`  
All: `slantRangeOffset`, `slantRangeOffsetVariance`, `alongTrackOffset`,
`alongTrackOffsetVariance`, `crossOffsetVariance`, `correlationSurfacePeak`,
`snr`, `digitalElevationModel`

### Key attributes

| Attribute | Description |
|---|---|
| `OffsetRangeSize` | Range dimension of offset map |
| `OffsetAzimuthSize` | Azimuth dimension of offset map |
| `OffsetRangePixelSize` | Effective range pixel size (m) |
| `OffsetAzimuthPixelSize` | Effective azimuth pixel size (m) |
| `r0`, `a0` | First range/azimuth sample offsets |
| `deltaR`, `deltaA` | Range/azimuth decimation factors |
| `layers` | List of loaded layer numbers |

### Writing offsets

`writeData` creates per-layer files with suffixes:

| Field | Suffix |
|---|---|
| `slantRangeOffset` | `.dr` |
| `alongTrackOffset` | `.da` |
| `slantRangeOffsetVariance` | `.vr` |
| `alongTrackOffsetVariance` | `.va` |
| `crossOffsetVariance` | `.vc` |
| `correlationSurfacePeak` | `.cc` |
| `snr` | `.snr` |

Each layer gets its own VRT: `filenameRoot.layer1.vrt`, `filenameRoot.layer2.vrt`, …

Pass `scaleToPixels=True` to convert metric offsets to pixel units.

---

## Common Range-Doppler attributes (all four classes)

These are available on every RD product after `openHDF`:

| Attribute | Description |
|---|---|
| `referenceOrbit` | Reference orbit number |
| `track` | Track number |
| `frame` | Frame number |
| `datetime` | `datetime` object for scene centre time |
| `centerLat`, `centerLon` | Scene centre lat/lon (degrees) |
| `LookDirection` | `'left'` or `'right'` |
| `PassType` | `'ascending'` or `'descending'` |
| `Wavelength` | Radar wavelength (m) |
| `PRF` | Effective PRF (Hz) |
| `orbit` | `nisarOrbit` object with state vectors |
| `polarization` | Active polarization string |
| `polarizations` | All available polarizations |
| `epsg` | EPSG code (4326 for lat/lon reference) |

---

## Geotransform convention

`getGeoTransform(tiff=False)` returns `[rOrigin, dR, 0, aOrigin, 0, dA]` where:
- `rOrigin` = near-range slant range − dR/2
- `dR` = range pixel spacing (positive, increases east/away)
- `aOrigin` = first zero-Doppler time − dA/2
- `dA` = zero-Doppler time spacing (**positive**, azimuth line 0 = earliest time)

With `tiff=True` the origin shifts to the last azimuth time and `dA` is
negated, giving a north-up GDAL convention.
