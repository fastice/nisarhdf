# showimage

Interactive viewer for NISAR HDF5 products, VRT, and GeoTIFF images. Displays 1–3 same-size images side by side with a floating control palette for colormap selection, modulo display, and per-pane band switching. Band data is loaded lazily — only read when first displayed.

## Entry points

| Command | Equivalent to |
|---|---|
| `showimage FILE [FILE ...]` | General viewer |
| `showvel FILE` | `showimage FILE --vel` — displays speed from a 2-band vx/vy file |
| `showoffsets FILE` | `showimage FILE --bands RangeOffsets AzimuthOffsets Correlation` |

## Usage

```
showimage [options] FILE [FILE ...]
```

**Positional:**

| Argument | Description |
|---|---|
| `FILE` | Input file(s) — up to 3; `.h5`/`.he5`/`.hdf5` for NISAR, `.vrt`/`.tif`/`.tiff` for GDAL |

**Display options:**

| Flag | Default | Description |
|---|---|---|
| `--cmap CMAP` | `gray` | Colormap name (matplotlib) |
| `--vmin X` | 2nd percentile | Lower clip value |
| `--vmax X` | 98th percentile | Upper clip value |
| `--mod X` | off | Display image modulo X (useful for phase wrapping) |
| `--decFactor N` | auto-fit screen | Decimation factor |
| `--fullRes` | — | Display at full resolution (equivalent to `--decFactor 1`) |

**Velocity display (GDAL only):**

| Flag | Description |
|---|---|
| `--vel` | Read bands 1+2 as vx, vy; display speed = √(vx²+vy²) mod 100 |
| `--log` | With `--vel`: log-scaled HSV speed rendering (1–3000 m/yr) |

**Band selection (GDAL only):**

| Flag | Description |
|---|---|
| `--bands BAND [BAND ...]` | Show 1–3 named bands from a single multi-band file |

**NISAR HDF5 options:**

| Flag | Default | Description |
|---|---|---|
| `--freq FREQ` | `frequencyA` | Frequency band to display |
| `--pol POL` | auto-detect | Polarization (e.g. `HH`, `VV`) |
| `--noCache` | — | Disable decimated-band cache; reduces memory use at the cost of re-reading from disk on each band switch |

## NISAR HDF5 support

Pass one or more `.h5`/`.he5`/`.hdf5` NISAR product files directly. The viewer auto-detects the product type from the HDF5 path and lists all displayable fields as band-switch buttons in the palette.

### Supported products and fields

| Product | Fields displayed |
|---|---|
| **RIFG** | `wrappedInterferogram` (as phase), `coherenceMagnitude` |
| **RUNW** | `unwrappedPhase`, `coherenceMagnitude`, `connectedComponents`, `ionospherePhaseScreen`, `ionospherePhaseScreenUncertainty` |
| **ROFF** | per layer: `slantRangeOffset`, `alongTrackOffset`, `correlationSurfacePeak`, `snr`, `slantRangeOffsetVariance`, `alongTrackOffsetVariance`, `crossOffsetVariance` |
| **GUNW** | `unwrappedPhase`, `coherenceMagnitude`, `connectedComponents`, `ionospherePhaseScreen`, `ionospherePhaseScreenUncertainty` |
| **GOFF** | same per-layer fields as ROFF |
| **GCOV** | covariance terms (e.g. `HHHH`, `VVVV`) + `mask`, `numberOfLooks` |

Complex-valued fields (e.g. `wrappedInterferogram`) are displayed as phase (angle).

### Band caching

The decimated array for each band is cached in memory after the first display so that switching back to a previously viewed band is instant. Use `--noCache` to disable this on memory-constrained machines.

## Interactive controls

The floating palette (left window) provides:

- **Colormap selector** — live colormap switching
- **Modulo input** — enter a value and press Return to apply modulo display
- **Band buttons** — one button per available band/field; click to switch the pane
- **Pick mode** — click a pixel to read its value
- **Profile mode** — click two points to extract and plot a line profile
- **Col Plot / Row Plot** — extract a column or row profile

Multiple image panes scroll in sync.

## Examples

```bash
# Display a RIFG interferogram (shows wrapped phase; coherence button in palette)
showimage NISAR_L1_PR_RIFG_*.h5

# Display with a custom colormap and phase limits
showimage NISAR_L1_PR_RIFG_*.h5 --cmap seismic --vmin -3.14 --vmax 3.14

# Compare two NISAR products (must be same size)
showimage NISAR_L1_PR_RIFG_*.h5 NISAR_L1_PR_RUNW_*.h5

# Display a GOFF geocoded offset product, frequencyB
showimage NISAR_L2_PR_GOFF_*.h5 --freq frequencyB

# Velocity file (VRT with vx, vy bands)
showvel velocity.vrt

# Offset correlation quicklook
showoffsets offsets.vrt
```
