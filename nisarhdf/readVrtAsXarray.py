#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 09:02:52 2024

@author: ian
"""
import rioxarray


def readVrtAsXarray(filename, mask_and_scale=True):
    '''
    Read a vrt and save each band as named variabile. Note this requires the
    band name be save and metadata field "Description" for each band as is
    for vrts saved with nisarHDF functions

    Parameters
    ----------
    filename : str
        Vrt file.

    Returns
    -------
    xarray : rioxarray.
        Data as a rioxarray.
    '''
    xarray = rioxarray.open_rasterio(filename,
                                     band_as_variable=True,
                                     mask_and_scale=True)
    # extract band names from 'Description' in each band
   # bandNames = [getattr(xarray[b], 'Description') for b in xarray.data_vars]
    bandNames = []
    for b in xarray.data_vars:
        attrs = xarray[b].attrs
        name = attrs.get('Description', attrs.get('long_name', b))
        bandNames.append(name)
    # Rename the bands
    if len(set(bandNames)) == len(bandNames):
        xarray = xarray.rename(dict(zip(xarray.data_vars, bandNames)))
    return xarray
