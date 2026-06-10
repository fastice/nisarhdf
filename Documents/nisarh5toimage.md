# nisarh5toimage — Command Line Reference

Converts NISAR HDF5 products to GeoTIFF, COG, flat binary, or VRT-only output.
Installed as an entry point by `pip install`.

---

## Synopsis

```
nisarh5toimage  productName  [output]  [options]
```

| Argument | Description |
|---|---|
| `productName` | NISAR HDF5 file (local path, `s3://…`, or `https://…`) |
| `output` | Root name for output files (omit with `--info`) |

---

## Options

### Output format

| Flag | Default | Description |
|---|---|---|
| `--outputFormat {COG,GTiff,binary}` | `COG` | Cloud-optimised GeoTIFF, standard GeoTIFF, or flat binary |
| `--vrtOnly` | off | Write a `.vrt` pointing directly into the HDF5; **no data is extracted** |
| `--MSB` | off | Write binary files big-endian (MSB); default is LSB |
| `--noSuffix` | off | Omit band name from filename (single-band output only) |
| `--scale SCALE` | 1.0 | Multiply every band by this factor before writing |

### Product selection

| Flag | Default | Description |
|---|---|---|
| `--productFamily` | parsed from filename | Force product type: `RSLC`, `ROFF`, `RIFG`, `RUNW`, `GCOV`, `GUNW`, `GOFF` |
| `--frequencyB` | off | Use frequency B band instead of frequency A |
| `--polarization {HH,VV,HV,VH}` | first like-pol | Polarization for non-GCOV products |
| `--fields FIELD …` | product defaults | Fields to extract; use `all` for everything |
| `--layers {1,2,3} …` | `1 2 3` | Offset layers (ROFF/GOFF only) |
| `--wrapped` | off | Select wrapped interferogram (GUNW only) |

### GCOV-specific

| Flag | Default | Description |
|---|---|---|
| `--dB` | off | Convert power to dB (10 log₁₀) |
| `--sigma0` | off | Output sigma0 instead of gamma0 |

### ROFF/GOFF-specific

| Flag | Default | Description |
|---|---|---|
| `--scaleToPixels` | off | Convert metric offsets to pixel units |

### Resolution

| Flag | Default | Description |
|---|---|---|
| `--downsampleFactor N` | 1 | Reduce resolution by factor N (integer > 0) |
| `--quickLook` | off | Write PNG quick-look (uses default downsample for product) |

### Access

| Flag | Default | Description |
|---|---|---|
| `--ros3` | off | Use HDF5 ROS3 driver for S3/HTTPS (slower, less memory than buffered read) |
| `--conserveMem` | off | Stream from HDF5 without converting to numpy (slow, low memory) |

### Info

| Flag | Description |
|---|---|
| `--info` | Print metadata and available fields; no output written |

---

## Output files

### Binary / tiff output

For non-offset products, each field is written as a separate file:

```
output.unwrappedPhase           (binary)
output.unwrappedPhase.tif       (GeoTIFF/COG)
output.vrt                      (multi-band VRT linking all fields)
```

For offset products, files are per-layer with GrIMP suffixes:

```
output.layer1.dr                slant-range offset
output.layer1.da                along-track offset
output.layer1.vr                range offset variance
…
output.layer1.vrt               multi-band VRT for layer 1
output.layer2.vrt               multi-band VRT for layer 2
```

### VRT-only output (`--vrtOnly`)

Only a `.vrt` is written; the HDF5 file is not opened for data.
The VRT contains `SimpleSource` elements pointing to the HDF5 bands
using the GDAL HDF5 subdataset path format:

```
HDF5:"/abs/path/to/file.h5"://science/LSAR/RSLC/swaths/frequencyA/HH
```

No binary extraction, no copy — the original HDF5 is the data store.

---

## Examples

### Inspect a product

```bash
nisarh5toimage  NISAR_L1_PR_RUNW_…h5  --info
```

### RUNW → COG (default)

```bash
nisarh5toimage  NISAR_L1_PR_RUNW_…h5  output_root
# writes: output_root.unwrappedPhase.tif  output_root.vrt
```

### RUNW → binary, all fields

```bash
nisarh5toimage  NISAR_L1_PR_RUNW_…h5  output_root \
    --outputFormat binary  --fields all
```

### RSLC → VRT-only (no extraction)

```bash
nisarh5toimage  NISAR_L1_PR_RSLC_…h5  output_root  --vrtOnly
# writes: output_root.vrt  (points to HH band in the HDF5)
# Verify: gdalinfo output_root.vrt
```

### GCOV → downsampled COG in dB

```bash
nisarh5toimage  NISAR_L2_PR_GCOV_…h5  output_root \
    --fields HHHH VVVV  --dB  --downsampleFactor 8
```

### GCOV quick-look PNG

```bash
nisarh5toimage  NISAR_L2_PR_GCOV_…h5  output_root  --quickLook
```

### ROFF layer 1, scale to pixels, binary

```bash
nisarh5toimage  NISAR_L1_PR_ROFF_…h5  output_root \
    --outputFormat binary  --scaleToPixels  --layers 1
```

### GUNW wrapped interferogram

```bash
nisarh5toimage  NISAR_L2_PR_GUNW_…h5  output_root  --wrapped
```

### Remote S3 file (low memory)

```bash
nisarh5toimage  s3://nisar-data/…/NISAR_L1_PR_RUNW_…h5  output_root \
    --ros3
```

---

## VRT structure

The companion `.vrt` includes:

- `<GeoTransform>` — range origin / pixel size for RD products; x/y origin for GEO
- `<Metadata>` — all parsed product parameters (orbit, timing, geometry, state vectors)
- One `<VRTRasterBand>` per field

For **binary output** bands use `VRTRawRasterBand` with `relativeToVRT="1"`.  
For **VRT-only** bands use `VRTRasterBand` + `SimpleSource` with the absolute
HDF5 subdataset path.

Reading the VRT back into Python:

```python
import nisarhdf
data = nisarhdf.readVrtAsXarray('output_root.vrt')
phase = data['unwrappedPhase']   # rioxarray DataArray
```

---

## Exit codes

| Code | Meaning |
|---|---|
| 0 | Success |
| 1 | Fatal error (bad arguments, missing file, etc.) |
