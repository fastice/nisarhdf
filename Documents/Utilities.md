# Utility Modules

---

## writeMultiBandVrt

**Source:** `writeMultiBandVrt.py`  
**Exported as:** `nisarhdf.writeMultiBandVrt`

Writes a GDAL VRT that references one or more source files as separate bands.
Used internally by `writeData` / `_writeNonOffsetData` after writing binary or
tiff files.

```python
from nisarhdf import writeMultiBandVrt

writeMultiBandVrt(
    newVRTFile,            # output .vrt path
    xSize, ySize,          # raster dimensions
    sourceFiles,           # list of source filenames
    descriptions,          # list of band description strings
    byteOrder='MSB',       # 'MSB' or 'LSB' (binary sources)
    eType=gdal.GDT_Float32,  # GDAL type or list of types
    geoTransform=[…],      # standard 6-element geotransform
    metaData={},           # dict of dataset-level metadata
    epsg=None,             # EPSG code for SRS
    noDataValue=-2.0e9,    # nodata value or list
    tiff=True,             # True → _createTiffVrt, False → _createBinaryVrt
    scales=None,           # optional per-band scale factors
    offsets=None           # optional per-band offsets
)
```

**Tiff sources** use `gdal.BuildVRT` (separate=True).  
**Binary sources** use the GDAL VRT driver directly with `VRTRawRasterBand`.

---

## writeVrtOnly

**Source:** `writeVrtOnly.py`

Writes a VRT that points **directly at HDF5 bands**, without extracting data.
Used by `writeData(vrtOnly=True)` and `nisarh5toimage --vrtOnly`.

### Public functions

```python
from nisarhdf.writeVrtOnly import getHDF5BandPath, writeVrtOnly
```

#### `getHDF5BandPath(productObj, field, layer=None) → str`

Returns the full internal HDF5 path for a field, e.g.:

```
/science/LSAR/RSLC/swaths/frequencyA/HH
/science/LSAR/RUNW/swaths/frequencyA/interferogram/HH/unwrappedPhase
/science/LSAR/ROFF/swaths/frequencyA/pixelOffsets/HH/layer1/slantRangeOffset
```

| Product | Path pattern |
|---|---|
| RSLC / GCOV / GSLC | `/science/{sar}/{product}/{bands}/{frequency}/{field}` |
| Others, no layer | `…/{productType}/{polarization}/{field}` |
| ROFF/GOFF with layer | `…/{polarization}/{layer}/{field}` |

#### `writeVrtOnly(productObj, filenameRoot, fields, vrtFile=None, layer=None)`

Writes a `.vrt` whose bands are `VRTRasterBand` + `SimpleSource` pointing to
GDAL HDF5 subdataset strings:

```
HDF5:"/abs/path/file.h5"://science/LSAR/RSLC/swaths/frequencyA/HH
```

Geotransform convention:
- RD products (`coordType='RD'`): `tiff=False` → positive dy (GrIMP binary)
- GEO products (`coordType='GEO'`): `tiff=True` → north-up convention

**Note:** remote (s3:// / https://) paths are embedded as-is; GDAL HDF5
virtual-I/O support is required for them to be read back.

---

## readVrtAsXarray

**Source:** `readVrtAsXarray.py`  
**Exported as:** `nisarhdf.readVrtAsXarray`

Reads a nisarhdf-produced `.vrt` back into an `xarray` / `rioxarray` dataset
with bands named by their `Description` metadata item.

```python
from nisarhdf import readVrtAsXarray

ds = readVrtAsXarray('output_root.vrt', mask_and_scale=True)
phase = ds['unwrappedPhase']    # rioxarray DataArray
coh   = ds['coherenceMagnitude']
```

Requires the VRT to have `Description` set per band (all nisarhdf-written VRTs do).

---

## nisarOrbit

**Source:** `nisarOrbit.py`

Stores and interpolates orbital state vectors.  Instantiated automatically
during `openHDF`; available as `productObj.orbit`.

```python
orbit = productObj.orbit

# State vector time series
orbit.position          # np.ndarray [N, 3]  ECEF position (m)
orbit.velocity          # np.ndarray [N, 3]  ECEF velocity (m/s)
orbit.TimeOfFirstStateVector   # seconds of day
orbit.StateVectorInterval      # spacing between state vectors (s)
orbit.NumberOfStateVectors

# Interpolated position / velocity at any time
pos, vel = orbit.interpolate(zdTime)   # zdTime in seconds of day
```

State vectors can be overridden from an external ESA/XML orbit file:

```python
myRUNW.openHDF('product.h5',
               referenceOrbitXML='S1A_OPER_AUX_POEORB_…xml')
```

---

## multilookRSLC

**Source:** `multilookRSLC.py`  
**Entry point:** `multilookRSLC`

Command-line tool to multilook an RSLC to a power image and write it to
COG, GeoTIFF, or binary.

```
multilookRSLC  RSLCFile  output  [options]

  --frequencyB            Use frequency B
  --dB                    Output in dB
  --polarization POL      HH (default) | VV | HV | VH
  --outputFormat FORMAT   COG | GTiff | binary
  --downsampleFactor N    Average N × N pixels
```

---

## s3Listing

**Source:** `s3Listing.py`  
**Entry point:** `s3Listing`

Recursive S3 bucket listing with filtering, date range, and tree/flat display.

```bash
s3Listing  s3://nisar-data/NISAR/   \
    --fileEndsWith .h5              \
    --nameFilter "RUNW*"            \
    --createdAfter 2025-01-01       \
    --flat
```

| Flag | Description |
|---|---|
| `--nameFilter PATTERN` | Substring or glob (e.g. `"7700*SH"`) |
| `--excludeName TEXT` | Exclude paths containing this text |
| `--fileEndsWith EXT` | File extension filter (e.g. `.h5`) |
| `--createdAfter DATE` | ISO date `YYYY-MM-DD` |
| `--createdBefore DATE` | ISO date `YYYY-MM-DD` |
| `--nonRecursive` | List top-level keys only |
| `--flat` | Print flat list instead of tree |
| `--long` | Include size and timestamp |
| `--profile PROFILE` | AWS profile name |

---

## formatGeojson

**Source:** `formatGeojson.py`  
**Exported as:** `nisarhdf.formatGeojson`

Formats a scene footprint GeoJSON feature for GrIMP geodat files.

```python
from nisarhdf import formatGeojson
geojson_str = formatGeojson(corners, properties_dict)
```

---

## Plotting tools

**Source:** `nisarhdfPlottingTools.py`  
**Exported as:** `nisarhdf.autoScaleRange`, `nisarhdf.colorBar`, `nisarhdf.createDivider`

Small matplotlib helpers used in display notebooks.

```python
from nisarhdf import autoScaleRange, colorBar, createDivider

vmin, vmax = autoScaleRange(data, percentile=98)

colorBar(pos, ax, label='Phase (rad)',
         colorBarPosition='right',
         colorBarSize='5%',
         colorBarPad=0.05)

# Create a divider for adding a colorbar axis
div = createDivider(ax, colorBarPosition='right', colorBarSize='5%')
```
