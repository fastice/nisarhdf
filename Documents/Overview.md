# nisarhdf — Package Overview

Python package for reading, processing, and reformatting NISAR HDF5 products.
Designed primarily for the **GrIMP** (Greenland Ice Mapping Project) workflow but
the core classes are product-agnostic and useful for any NISAR data work.

---

## Installation

```bash
cd /path/to/nisarhdf
pip install -e .
```

Dependencies: `h5py`, `numpy`, `scipy`, `gdal` (osgeo), `rioxarray`, `pyproj`,
`boto3` (for S3 access).

---

## What the package does

- Opens NISAR HDF5 files from local disk, S3, or HTTPS and exposes bands as
  plain numpy arrays (e.g. `myRUNW.unwrappedPhase`).
- Parses all essential metadata into Python attributes (orbit, frame, timing,
  pixel spacing, incidence angle, wavelength, …).
- Interpolates HDF5 metadata *cubes* (incidenceAngle, baseline, los/along-track
  unit vectors, ground-track velocity) onto the data grid.
- Applies masks and ionospheric corrections.
- Writes bands to flat binary, GeoTIFF, COG, or PNG with companion `.vrt`.
- Alternatively writes a lightweight `.vrt` that points directly into the HDF5
  (no data extraction — see `--vrtOnly` / `writeVrtOnly`).
- Provides the `nisarh5toimage` command-line tool for batch conversion.

---

## Supported products

| Class | Product | Coordinates | lookType |
|---|---|---|---|
| `nisarRSLCHDF` | RSLC | Range-Doppler | `SLC` |
| `nisarRUNWHDF` | RUNW | Range-Doppler | `ML` |
| `nisarRIFGHDF` | RIFG | Range-Doppler | `ML` |
| `nisarROFFHDF` | ROFF | Range-Doppler | `Offset` |
| `nisarGCOVHDF` | GCOV | Geocoded (x/y) | `None` |
| `nisarGUNWHDF` | GUNW | Geocoded (x/y) | `ML` |
| `nisarGOFFHDF` | GOFF | Geocoded (x/y) | `Offset` |

---

## Class hierarchy

```
nisarBaseHDF  (nisarBaseHDF.py)
│
├── nisarBaseRangeDopplerHDF  (nisarBaseRangeDopplerHDF.py)
│   │   coordType = 'RD'  — coordinates are slant range & zero-Doppler time
│   ├── nisarRSLCHDF
│   ├── nisarRUNWHDF
│   ├── nisarRIFGHDF
│   └── nisarROFFHDF
│
└── nisarBaseGeocodedHDF  (nisarBaseGeocodedHDF.py)
    │   coordType = 'GEO'  — coordinates are projected x/y (polar stereo or UTM)
    ├── nisarGCOVHDF
    ├── nisarGUNWHDF
    └── nisarGOFFHDF
```

---

## Quick start

```python
import nisarhdf

# --- RUNW (unwrapped phase) ---
myRUNW = nisarhdf.nisarRUNWHDF()
myRUNW.openHDF('NISAR_L1_PR_RUNW_…h5')

print(myRUNW.NumberRangeLooks, myRUNW.PRF)
phase = myRUNW.unwrappedPhase          # numpy float32 array
coh   = myRUNW.coherenceMagnitude

# Apply coherence mask and ionosphere correction
myRUNW.maskPhase()
myRUNW.ionosphereCorrectPhase()

# Save to binary + VRT
myRUNW.writeData('output_root', outputFormat='binary')

# --- RSLC ---
myRSLC = nisarhdf.nisarRSLCHDF()
myRSLC.openHDF('NISAR_L1_PR_RSLC_…h5', fields=['HH'])
slc = myRSLC.HH    # complex64 numpy array

# VRT-only (no data extraction)
myRSLC2 = nisarhdf.nisarRSLCHDF()
myRSLC2.openHDF('NISAR_L1_PR_RSLC_…h5', noLoadData=True, fields=['HH'])
myRSLC2.writeData('output_root', vrtOnly=True)
```

---

## Coordinate conventions

**Range-Doppler products (RSLC, RUNW, RIFG, ROFF)**

The geotransform uses:
- x-axis: slant range (metres), positive east  
- y-axis: zero-Doppler time (seconds of day), **positive downward** (dy > 0)

This is the GrIMP binary convention (`tiff=False`). When writing to GeoTIFF the
y-axis is flipped (dy < 0) to match the GDAL/GIS north-up expectation.

**Geocoded products (GCOV, GUNW, GOFF)**

Standard projected x/y coordinates (polar stereographic or UTM depending on
hemisphere). EPSG code is read from the HDF5 and stored as `self.epsg`.
dy follows the tiff convention (negative = north-up) for both GeoTIFF and VRT
outputs.

---

## Further reading

| Topic | Document |
|---|---|
| Range-Doppler product classes (RSLC, RUNW, RIFG, ROFF) | [RangeDopplerProducts.md](RangeDopplerProducts.md) |
| Geocoded product classes (GCOV, GUNW, GOFF) | [GeocodedProducts.md](GeocodedProducts.md) |
| Base class method reference | [BaseClasses.md](BaseClasses.md) |
| `nisarh5toimage` CLI tool | [nisarh5toimage.md](nisarh5toimage.md) |
| Utility modules | [Utilities.md](Utilities.md) |
