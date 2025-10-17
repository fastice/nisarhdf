#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 08:25:27 2024

@author: ian
"""
from osgeo import gdal, osr
import os
import numpy as np

# gdalTypes = {'Float32': gdal.GDT_Float32,
#              'Float64': gdal.GDT_Float64,
#              'Complex64':  gdal.GDT_CFloat32,
#              'UInt16': gdal.GDT_UInt16,
#              'Int8': gdal.GDT_Byte,
#              'UInt8': gdal.GDT_Byte}


def writeMultiBandVrt(newVRTFile, xSize, ySize, sourceFiles, descriptions,
                      byteOrder='MSB', eType=gdal.GDT_Float32,
                      geoTransform=[-0.5, 1., 0., -0.5, 0., 1.], metaData={},
                      epsg=None, noDataValue=-2.0e9, tiff=True):
    '''
    Write a vrt for for several products as bands.
    Parameters
    ----------
    newVRTFile : str
        Name for vrt file being produced.
    sourceFiles : list
        List of files to serve as bands in vrt.
    descriptions : list
        List of descriptions for the bands.
    byteOrder : str, optional
        Byte order (MSB, LSB). If none, the default is use value from the
        metaData, if none there default to MSMB. The default is None.
    eType : GDT, optional
        Gdal data type or list of data types. The default is gdal.GDT_Float32.
    geoTransform : list, optional
        The geoTransform. The default is [-0.5, 1., 0., -0.5, 0., 1.] for
        array with no coordinates.
    metaData : dict, optional
        Dictionary with meta data values. The default is None.
    epsg : str, optional
        Seet to epsg value to set SRS. The default is None.
    noDataValue : int, float as appropriate, optional
        The nodata value. The default is -2.0e9.

    Returns
    -------
    None.

    '''
    if not isinstance(eType, list):
        eType = [eType] * len(sourceFiles)
    if not isinstance(noDataValue, list):
        noDataValue = [noDataValue] * len(sourceFiles)
    if type(sourceFiles) is not list:
        sourceFiles = [sourceFiles]
    if type(descriptions) is not list:
        descriptions = [descriptions]
    if byteOrder not in ['MSB', 'LSB']:
        print(f"Invalid byte order {byteOrder}")
        return
    # Kill any old file
    if os.path.exists(newVRTFile):
        os.remove(newVRTFile)
    # Create VRT
    if tiff:
        _createTiffVrt(newVRTFile, sourceFiles, descriptions, noDataValue,
                       geoTransform=geoTransform,
                       metaData=metaData,
                       epsg=epsg)
    else:
        _createBinaryVrt(newVRTFile, xSize, ySize, sourceFiles, descriptions,
                         eType, noDataValue,
                         byteOrder=byteOrder,
                         geoTransform=geoTransform,
                         metaData=metaData,
                         epsg=epsg)


def _createTiffVrt(newVRTFile, sourceFiles, descriptions, noDataValues,
                   geoTransform=[-0.5, 1., 0., -0.5, 0., 1.], metaData=None,
                   epsg=None):
    '''
    Write a tiff vrt using gdal.BuildVRT

    Parameters
    ----------
    newVRTFile : str
        Vrt files.
    sourceFiles : list of str
        Source file names.
    descriptions : list of str
        Description for each band.
    noDataValues :  list
        No data values for each band.
    geoTransform : list of floats, optional
        Standard geotransform. The default is [-0.5, 1., 0., -0.5, 0., 1.].
    metaData : dict, optional
        Dict with meta data that applies to all bands. The default is None.
    epsg : int, optional
        EPSG Code. The default is None.

    Returns
    -------
    None.

    '''
    # Build the vrt
    try:
        # This should work, but fails with gdal.3.9.1. It seems to default to
        # relativeToVRT=1, so maybe ok to leave out, but explicitly set if 
        # possible.
        options = gdal.BuildVRTOptions(options=["relativeToVRT=1"], separate=True)
    except Exception:
        options = gdal.BuildVRTOptions(separate=True)
    ds = gdal.BuildVRT(newVRTFile, sourceFiles, options=options)
    # print('meta', metaData)
    if bool(metaData):
        ds.SetMetadata(metaData)
    # This should be unnecesary since tiff should have info.
    if epsg is not None:
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(epsg)
        ds.SetSpatialRef(sr)
    #
    nBands = len(sourceFiles)
    # Update bands
    for description, bandNumber, noDataValue in \
            zip(descriptions, range(1, nBands + 1), noDataValues):
        #
        band = ds.GetRasterBand(bandNumber)
        band.SetMetadataItem("Description", description)
        if not isinstance(noDataValue, np.complex64):
            band.SetNoDataValue(noDataValue)
    ds.FlushCache()
    ds = None


def _createBinaryVrt(newVRTFile, xSize, ySize, sourceFiles, descriptions,
                     eTypes, noDataValues,
                     byteOrder='LSB',
                     geoTransform=[-0.5, 1., 0., -0.5, 0., 1.],
                     metaData=None,
                     epsg=None):
    '''
    Write a tiff vrt using vrt driver directly since BuildVRT doesn't appear
    to work with binary files.

    Parameters
    ----------
    newVRTFile : str
        Vrt files.
    sourceFiles : list of str
        Source file names.
    descriptions : list of str
        Description for each band.
    eTypes : list
        Data types in gdal codes
    noDataValues :  list
        No data values for each band.
    byteOrder : str
        Byte order. The default is 'LSB'
    geoTransform : list of floats, optional
        Standard geotransform. The default is [-0.5, 1., 0., -0.5, 0., 1.].
    metaData : dict, optional
        Dict with meta data that applies to all bands. The default is None.
    epsg : int, optional
        EPSG Code. The default is None.

    Returns
      -------
      None.
    '''
    bands = len(sourceFiles)
    drv = gdal.GetDriverByName("VRT")
    vrt = drv.Create(newVRTFile, xSize, ySize, bands=0)
    vrt.SetGeoTransform(geoTransform)
    #
    if epsg is not None:
        sr = osr.SpatialReference()
        sr.ImportFromEPSG(epsg)
        vrt.SetSpatialRef(sr)
    if metaData is not None:
        vrt.SetMetadata(metaData)
    # Loop to add bands
    for sourceFile, description, bandNumber, dataType, noDataValue in \
            zip(sourceFiles, descriptions, range(1, bands + 1), eTypes,
                noDataValues):
        # setup options with filename
        options = [f"SourceFilename={os.path.basename(sourceFile)}",
                   "relativeToVRT=1",
                   "subclass=VRTRawRasterBand"]
        if dataType not in [gdal.GDT_Byte]:
            options.append(f"BYTEORDER={byteOrder}")
        options.append(bytes(0))
        #
        vrt.AddBand(dataType, options=options)
        band = vrt.GetRasterBand(bandNumber)
        band.SetMetadataItem("Description", description)
        if not isinstance(noDataValue, np.complex64):
            band.SetNoDataValue(noDataValue)
    # Close the vrt
    vrt = None