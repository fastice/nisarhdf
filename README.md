# nisarHDF

This repository contains python code to read the NISAR HDF formats and parse out the necessary parameters and reformatting necessary for the Greenland Ice Mapping Project (*GrIMP*) workflow. It also provide a command line utility (`nisarh5toimage`) for reformatting and extracting NISAR HDF5 data. Most of the functionality though is independent of *GrIMP* and may prove useful for anyone working with NISAR HDF data.

## Python Package
In particular, it simplifies access to NISAR data by allowing:
  - Access data (e.g., `myRUNW.unwrappedPhase`) as numpy variables and other commonly used parameters (e.g., `myGUNW.NumberRangeLooks`) rather than hdf structures,
  - Mask data (e.g., remove  data with connected component values of zero),
  - Apply ionospheric corrections to phase data,
  - Interpolate metadata cubes (e.g., incidenceAngle and baseline), including mapping these variables to the same grid as the data products.
  - Save data and interpolated metadata cubes to binary or geotiff files, and
  - Write a .vrt file to link multiple data layers so they can be read back in as `xarray`.

At present classes have been developed for **RSLC**, **ROFF**, **RUNW**, **RIFG**, **GOFF**, **GUNW**, and **GCOV** products. A tutorial notebook for each of these data types is included in the [Notebook](https://github.com/fastice/nisarhdf/blob/main/Notebooks) folder of this repository.

## Documentation

- [Overview](Documents/Overview.md) — package overview and class hierarchy
- [BaseClasses](Documents/BaseClasses.md) — base HDF reader classes
- [RangeDopplerProducts](Documents/RangeDopplerProducts.md) — RSLC, ROFF, RUNW, RIFG classes
- [GeocodedProducts](Documents/GeocodedProducts.md) — GOFF, GUNW, GCOV classes
- [nisarh5toimage](Documents/nisarh5toimage.md) — command-line HDF5 extraction utility
- [showimage](Documents/showimage.md) — interactive viewer for NISAR HDF5, VRT, and GeoTIFF images
- [Utilities](Documents/Utilities.md) — supporting utility functions

## Command line utility nisarh5toimage

 While HDF5 is a powerful format for archiving, other formats such as tiff and simple binary can be easier for many users to work with. To facilitate working with NISAR data for these users, the `nisarh5toimage` command line utility can:

- Provide a quick summary of what's in the file (bands etc) along with essential meta data (e.g., *time, orbit, frame, track, center lat/lon*, etc.)
- Pull the individual bands from NISAR hdfs and write them to tiff or binary files.
- Produce an accompanying `.vrt` field that makes it easy to import multiband products into Python (e.g., `rioxarray`) or GIS packages.
- Downsample products (useful for large hi-res GCOV products).
- Produce quickLook PNG files with `.vrt` for geospatial information.


The [nisarh5toimageTutorial](https://github.com/fastice/nisarhdf/blob/main/Notebooks/nisarh5toimageTutorial.ipynb) notebook in the [Notebook](https://github.com/fastice/nisarhdf/blob/main/Notebooks) directory of this respository includes several examples of how to use this program.
